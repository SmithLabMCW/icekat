import re
from decimal import Decimal
import math
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from uncertainties import ufloat

def logWithZeros(x):
    '''
    return log10 of array that may contain zeros
    '''
    out = []
    if len(x) > 0:
        for xi in x:
            if xi == 0.:
                out.append(0.)
            else:
                out.append(np.log10(xi))
    return np.array(out)

def linear(x, m, b):
    '''
    straight line
    '''
    return m*x + b

def logarithmic(x, yo, b, to):
    '''
    logarithmic equation from Lu & Fei et. al, 2003
    '''
    return yo + b*np.log(1 + x*to)

def mmfit(x, km, vmax):
    '''
    Michaelis Menten equation
    '''
    return vmax * x / (km + x)

def icfit(x, bottom, top, slope, p50):
    '''
    IC50 equation
    '''
    return bottom + (top-bottom)/(1+10**((-p50-x)*slope))

def spline_fit(x, y):
    
    x, y = x.values, y.values
    spline = UnivariateSpline(x, y)(x)
    derivative = np.abs(np.diff(spline)/np.diff(x))
    threshold = 0.7*(np.max(derivative) - np.min(derivative)) + np.min(derivative)
    try:
        indices = np.where(derivative > threshold)[0]
    except:
        indices = []
    while len(indices) < 4:
        threshold = threshold*0.9
        try:
            indices = np.where(derivative > threshold)[0]
        except:
            indices = []
    xi, yi = x[indices], y[indices]
    df = pd.DataFrame(data={'x' : xi, 'y' : yi}).sort_values('x').dropna()
    xi, yi = df.x, df.y
    popt, pcov = curve_fit(linear, xi, yi)
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(np.min(x), np.max(x), len(spline))
    yfit = linear(xfit, *popt)
    fit_dict = { 'x' : x, 'y' : y,
                    'rate' : np.abs(popt[0]), 'error' : perr[0],
                    'xfit' : xfit,
                    'yfit' : yfit, 'resi' : np.array(yfit) - y }

    return fit_dict

def linear_fit(x, y):

    popt, pcov = curve_fit(linear, x, y)
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(np.min(x), np.max(x), len(x))
    yfit = linear(xfit, *popt)
    fit_dict = { 'x' : x, 'y' : y,
                    'rate' : np.abs(popt[0]), 'error' : perr[0],
                    'xfit' : xfit,
                    'yfit' : yfit, 'resi' : np.array(yfit) - y }

    return fit_dict

def logarithmic_fit(x, y):
    
    popt, pcov = curve_fit(logarithmic, x, y, maxfev=100000)
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(np.min(x), np.max(x), len(x))
    yfit = logarithmic(xfit, *popt)
    yerr = logarithmic(xfit, *perr)
    fit_dict = { 'x' : x, 'y' : y,
                    'rate' : np.abs(np.array(np.diff(yfit)/np.diff(xfit))[0]),
                    'error' : np.abs(np.array(np.diff(yerr)/np.diff(xfit))[0]),
                    'xfit' : xfit,
                    'yfit' : yfit, 'resi' : np.array(yfit) - y }

    return fit_dict

class progress_curve(object):

    def __init__(self, dataframe):

        self.data = dataframe

    def spline(self, start, end):

        df = self.data[(self.data[self.data.columns[0]] >= float(start)) &
                       (self.data[self.data.columns[0]] <= float(end))].dropna(axis=0)
        x, y = df[df.columns[0]], df[df.columns[1]]
        spline = spline_fit(x, y)
        self.spline = spline
        
        return self.spline

    def linear(self, start, end):

        df = self.data[(self.data[self.data.columns[0]] >= float(start)) &
                       (self.data[self.data.columns[0]] <= float(end))].dropna(axis=0)
        x, y = df[df.columns[0]], df[df.columns[1]]
        linear = linear_fit(x, y)
        self.linear = linear
        
        return self.linear
    
    def logarithmic(self, start, end, offset):
        
        df = self.data[(self.data[self.data.columns[0]] >= float(start)) &
                       (self.data[self.data.columns[0]] <= float(end))].dropna(axis=0)
        x, y = df[df.columns[0]], df[df.columns[1]]
        try:
            x = x + offset
        except:
            pass
        logarithmic = logarithmic_fit(x, y)
        self.logarithmic = logarithmic

        return self.logarithmic

class kinetic_model(object):

    def __init__(self, dictionary):
        self.dict = dictionary

    def model(self, subtract, transform, threshold, bottom, top, slope, scalex):
        
        result = {}
        df = pd.DataFrame()
        for s in self.dict:
            if type(self.dict[s]) == progress_curve:
                if len(re.findall(r'-?\d+\.?\d*', str(s))) > 0:
                    x = np.float(re.findall(r'-?\d+\.?\d*', str(s))[0])
                else:
                    x = 0.0
                if self.dict[s+'_fit'] == 0:
                    sdf = self.dict[s].spline
                elif self.dict[s+'_fit'] == 1:
                    sdf = self.dict[s].linear
                else:
                    sdf = self.dict[s].logarithmic
                df.at[s, 'rate'] = sdf['rate']
                df.at[s, 'error'] = sdf['error']
                df.at[s, 'x'] = x
        df = df.sort_values(['x'])
        uRates = [ufloat(ur, ue) for ur, ue in zip(df['rate'], df['error'])]
        df['uRates'] = uRates
        try:
            df['uRates'] = df['uRates'] - df.loc[subtract]['uRates']
            df['rate'] = [ur.nominal_value for ur in df['uRates']]
            df['error'] = [ur.std_dev for ur in df['uRates']]
            df = df[df.index != subtract]
        except:
            pass
        try:
            x = np.array([ufloat(r, e) for r, e in zip(df['rate'], df['error'])])
            x = eval(transform)
            df['rate'] = [xi.nominal_value for xi in x]
            df['error'] = [xi.std_dev for xi in x]
        except:
            pass
        n, x, y, e = df.index.values, df['x'].values, df['rate'].values, df['error'].values
        result['n'] = n
        result['xt'] = x
        result['yt'] = ['%.2E' % Decimal(str(yi)) for yi in y]
        result['et'] = ['%.2E' % Decimal(str(ei)) for ei in e]
        result['l'] = y - e
        result['u'] = y + e
        if self.dict['model'] == 'Michaelis-Menten':
            include = np.where(x != 1e-23)[0]
            x, y, e = x[include], y[include], e[include]
            result['yp'] = y
            result['xp'] = x
            result['ep'] = e
            result['cp'] = ['grey']*len(result['xp'])
            result['ct'] = ['white']*len(result['xt'])
            try:
                xfit = np.linspace(np.min(x), np.max(x), 100)
                result['xfit'] = xfit
                popt_mm, pcov_mm = curve_fit(mmfit, x, y, sigma=e, absolute_sigma=True, maxfev=100000)
                perr_mm = np.sqrt(np.diag(pcov_mm))
                ymm = mmfit(xfit, *popt_mm)
                result['yfit'] = ymm
                result['Km'] = tuple(['%.2E' % Decimal(str(popt_mm[0])), 
                                        '%.2E' % Decimal(str(perr_mm[0]))])
                result['Vmax'] = tuple(['%.2E' % Decimal(str(popt_mm[1])), 
                                        '%.2E' % Decimal(str(perr_mm[1]))])
            except:
                result['xfit'] = []
                result['yfit'] = []
                result['Km'] = tuple([np.nan, np.nan])
                result['Vmax'] = tuple([np.nan, np.nan])
        elif self.dict['model'] == 'pEC50/pIC50':
            include = np.where(x != 1e-23)[0]
            x, y, e = x[include], y[include], e[include]
            include = np.where(x != 0.)[0]
            x, y, e = x[include], y[include], e[include]
            result['yp'] = y
            result['xp'] = x
            result['ep'] = e
            result['l'] = y - e
            result['u'] = y + e
            xn, yn, en = [], [], []
            for xin, yin, ein in zip(x, y, e):
                if xin != 0.:
                    xn.append(xin)
                    yn.append(yin)
                    en.append(ein)
            bounds = ([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf])
            for ix, b in enumerate([bottom, top, slope]):
                regex = re.findall(r'-?\d+\.?\d*', str(b))
                if len(regex) > 0:
                    regex = eval(b)
                    if regex > 0:
                        bounds[0][ix] = regex - regex/1e10
                        bounds[1][ix] = regex + regex/1e10
                    elif regex < 0:
                        bounds[0][ix] = regex + regex/1e10
                        bounds[1][ix] = regex - regex/1e10
                    elif regex == 0:
                        bounds[0][ix] = regex - 1e-10
                        bounds[1][ix] = regex + 1e-10
            result['cp'] = ['grey']*len(result['xp'])
            result['ct'] = ['white']*len(result['xt'])
            xfit = np.linspace(np.min(result['xp']), np.max(result['xp']), 100)
            try:
                popt_ic, pcov_ic = curve_fit(icfit, xn, yn, sigma=en, absolute_sigma=True, bounds=bounds, maxfev=100000)
                perr_ic = np.sqrt(np.abs(np.diag(pcov_ic)))
                for ix, b in enumerate([bottom, top, slope]):
                    regex = re.findall(r'-?\d+\.?\d*', str(b))
                    if len(regex) > 0:
                        popt_ic[ix] = eval(b)
                        perr_ic[ix] = 0.0
                yic = icfit(xfit, *popt_ic)
                result['xfit'] = xfit
                result['yfit'] = yic
                result['Bottom'] = np.array(['%.2E' % Decimal(str(popt_ic[0])),
                                                 '%.2E' % Decimal(str(perr_ic[0]))])
                result['Top'] = np.array(['%.2E' % Decimal(str(popt_ic[1])),
                                              '%.2E' % Decimal(str(perr_ic[1]))])
                result['Slope'] = np.array(['%.2E' % Decimal(str(popt_ic[2])), 
                                            '%.2E' % Decimal(str(perr_ic[2]))])
                result['p50'] = np.array(['%.2E' % Decimal(str(popt_ic[3])),
                                              '%.2E' % Decimal(str(perr_ic[3]))])
                result['cp'] = ['grey']*len(result['xp'])
                result['ct'] = ['white']*len(result['xt'])
            except:
                result['xfit'] = []
                result['yfit'] = []
                result['Bottom'] = [np.nan, np.nan]
                result['Top'] = [np.nan, np.nan]
                result['Slope'] = [np.nan, np.nan]
                result['p50'] = [np.nan, np.nan]               
            
        else:
            result['xt'] = np.linspace(1, len(n), len(n))
            result['yp'] = y
            result['xp'] = np.linspace(1, len(n), len(n))
            result['ep'] = e
            result['xfit'] = np.linspace(1, len(n), len(n))
            result['yfit'] = np.repeat(np.mean(result['yp']), len(n))
            std = np.std(result['yp'])
            avg = np.mean(result['yp'])
            color, colort = [], []
            for r in result['yp']:
                if r >= avg + std*threshold:
                    color.append('red')
                    colort.append('#EC7063')
                elif r <= avg - std*threshold:
                    color.append('blue')
                    colort.append('#5DADE2')
                else:
                    color.append('grey')
                    colort.append('white')
            result['cp'] = color
            result['ct'] = colort
        self.result = result
        if scalex == 1:
            result['xt'] = logWithZeros(result['xt'])
            result['xp'] = logWithZeros(result['xp'])
            result['xfit'] = np.linspace(np.min(result['xp']), np.max(result['xp']), 100)
            if self.dict['model'] == 'Michaelis-Menten':
                popt_mm, pcov_mm = curve_fit(mmfit, result['xp'], y, sigma=e, absolute_sigma=True, maxfev=100000)
                yic = mmfit(result['xfit'], *popt_mm)
                result['yfit'] = yic
            elif self.dict['model'] == 'pEC50/pIC50':
                popt_ic, pcov_ic = curve_fit(icfit, result['xp'], y, sigma=e, absolute_sigma=True, maxfev=100000)
                yic = icfit(result['xfit'], *popt_ic)
                result['yfit'] = yic
            
        return self.result
