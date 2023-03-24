import continuous_kinetics as ck
from os.path import dirname, join
from io import StringIO
import base64
import pandas as pd
import numpy as np
from bokeh import events
from bokeh.io import curdoc
from bokeh.layouts import row, column, layout
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, Div, BasicTickFormatter, Whisker, Range1d
from bokeh.models.widgets import RadioButtonGroup, Select, TextInput, Button, DataTable, TableColumn, RangeSlider, Slider, HTMLTemplateFormatter, CheckboxButtonGroup
from bokeh.plotting import figure

########## bokeh methods ##########
def widget_callback(attrname, old, new):

    start = float(start_time.value)
    stop = float(end_time.value)
    range_slider.value = (start, stop)

    update()

def sample_callback(attrname, old, new):

    start = float(range_slider.start)
    stop = float(range_slider.end)
    start_e = float(start_time.value)
    stop_e = float(end_time.value)

    if start_e != start or stop_e != stop:
        start_time.value = str(start)
        end_time.value = str(stop)

    update()

def slider_callback(attrname, old, new):

    start = float(range_slider.value[0])
    start_e = float(start_time.value)
    stop = float(range_slider.value[1])
    stop_e = float(end_time.value)
    step = float(range_slider.step)
    end = float(range_slider.end)

    if start >= stop -  4*step and stop <= end - 4*step:
        range_slider.value = (start, start + 4*step)

    elif start >= stop - 4*step and stop >= end - 4*step:
        range_slider.value = (stop - 4*step, stop)

    if start_e != start or stop_e != stop:
        start_time.value = str(start)
        end_time.value = str(stop)

    update()

def xbox_callback(attrname, old, new):

    start = float(range_slider.value[0])
    start_e = float(start_time.value)
    stop = float(range_slider.value[1])
    stop_e = float(end_time.value)
    begin = float(range_slider.start)
    step = float(range_slider.step)
    end = float(range_slider.end)

    if start_e != start or stop_e != stop:
        if start_e >= begin and stop_e <= end:
            range_slider.value = (start_e, stop_e)

    update()

def threshold_callback(attrname, old, new):

    if experiment_db['model'] == 'High-Throughput Screen':

        std = np.std(model_plot_source.data['yp'])
        avg = np.mean(model_plot_source.data['yp'])
        color, colort = [], []

        for r in model_plot_source.data['yp']:

            if r >= avg + std*threshold_slider.value:
                color.append('red')

            elif r <= avg - std*threshold_slider.value:
                color.append('blue')

            else:
                color.append('grey')

        for r in model_data_source.data['yt']:

            if float(r) >= avg + std*threshold_slider.value:
                colort.append('#EC7063')

            elif float(r) <= avg - std*threshold_slider.value:
                colort.append('#5DADE2')

            else:
                colort.append('white')

        model_plot_source.data['cp'] = color
        model_data_source.data['ct'] = colort

def update():

#    try:
        # get selections
        fit_routine = fit_button.active
        scalex = len(scalex_box.active)
        model_eq = model_select.value
        subtract = subtract_select.value
        sample = sample_select.value
        transform = transform_input.value
        offset = offset_input.value
        top = top_fix.value
        bottom = bottom_fix.value
        slope = slope_fix.value
        threshold = threshold_slider.value
        start = start_time.value
        end = end_time.value
        warning.visible = False
        warning_source.data = pd.DataFrame(data=dict(x=[], y=[], t=[]))
        circles.visible = False
        circles_source.data = pd.DataFrame(data=dict(x=[], y=[]))
        model.title.text = model_eq

        if 'x' in transform:
            raw.yaxis.axis_label = 'Concentration'
        else:
            raw.yaxis.axis_label = 'Concentration'

        # update database
        experiment_db[sample+'_fit'] = fit_routine
        experiment_db['model'] = model_eq
        pdf = experiment_df[[experiment_df.columns[0], sample]]
        experiment_db[sample] = ck.progress_curve(pdf, start, end)
        if fit_routine == 3:

            raw.title.text = "Schnell-Mendoza Fit"
            resi.title.text = "Schnell-Mendoza Fit Residuals"
            raw.title.text = "Schnell-Mendoza Fit"
            resi.title.text = "Schnell-Mendoza Fit Residuals"
            model_select.value = 'Michaelis-Menten'
            scalex_box.active = []

            model_data_source.data = pd.DataFrame(data = dict(
                                    xt = [], yt = [], et = [],
                                    n = [], ct = []
                                    )).to_dict('list')
            model_plot_source.data = pd.DataFrame(data = dict(xp = [],
                                    yp = [], u = [], l = [], ep = [],
                                    cp = [])).to_dict('list')
            ic_source.data = pd.DataFrame(data = dict(
                                    label = [], Bottom = [], Top = [],
                                    Slope = [], p50 = []
                                    )).to_dict('list')

            if 'x' in transform:

                raw_data, model_result, fit_data, varea_data = ck.sm_fit(experiment_db).fit(sample, transform, subtract)
                model_fit_source.data = pd.DataFrame(data = dict(
                                    x = model_result['xfit'],
                                    y = model_result['yfit']
                                    )).to_dict('list')
                varea_source.data = pd.DataFrame(data = dict(
                                    x = varea_data['x'],
                                    r1 = varea_data['r1'],
                                    r2 = varea_data['r2']
                                    )).to_dict('list')
                mm_source.data = pd.DataFrame(data = dict(
                                    label = ['Fit Value', 'Std. Error'],
                                    Km = fit_data['Km'],
                                    Vmax = fit_data['Vmax']
                                    ), index=['value', 'error']).to_dict('list')
                raw_source.data = pd.DataFrame(data = dict(
                                                     x = raw_data['x'], y = raw_data['y'],
                                                     yr = raw_data['resi'], yfit = raw_data['yfit'],
                                                    )).to_dict('list')

            else:
                warning.visible = True
                warning_source.data = pd.DataFrame(data=dict(x=[0], y=[0], t=['Please enter transform equation! \nMust convert signal to [substrate] \nin Schnell-Mendoza mode (e.g. via \nx/(ε *ℓ * [E]) for sample data). \nNote: this transform may need \nto be inverted through multiplying \nby -1 when analyzing experiments \nthat measure increasing product \nconcentration over time)']))
                circles.visible = True
                circles_source.data = pd.DataFrame(data=dict(x=[-.05, -.05, 1.6, 1.6], y=[0, 0.6, 0, 0.6]))
                raw.x_range = Range1d(-0.1, 2.5)
                raw_source.data = pd.DataFrame(data = dict(
                                                     x = [], y = [],
                                                     yr = [], yfit = [],
                                                    )).to_dict('list')
                model_fit_source.data = pd.DataFrame(data = dict(
                                    x = [],
                                    y = []
                                    )).to_dict('list')
                varea_source.data = pd.DataFrame(data = dict(
                                    x = [],
                                    r1 = [],
                                    r2 = []
                                    )).to_dict('list')
                mm_source.data = pd.DataFrame(data = dict(
                                    label = ['',''],
                                    Km = ['', ''],
                                    Vmax = ['', '']
                                    ), index=['value', 'error']).to_dict('list')

        else:

            raw.title.text = "Initial Rate Fit"
            resi.title.text = "Initial Rate Fit Residuals"
            varea_source.data = pd.DataFrame(data = dict( x = [], r1 = [], r2 = [] )).to_dict('list')
            raw.yaxis.axis_label = 'Signal'

            if fit_routine == 0:
                progress_data = experiment_db[sample].spline()
                experiment_db[sample].spline = progress_data

            elif fit_routine == 1:
                progress_data = experiment_db[sample].linear()
                experiment_db[sample].linear = progress_data

            elif fit_routine == 2:
                try:
                    offset = float(offset)
                except:
                        pass
                progress_data = experiment_db[sample].logarithmic(offset)
                experiment_db[sample].logarithmic = progress_data

            raw_source.data = pd.DataFrame(data = dict(
                                                    x = progress_data['x'],
                                                    y = progress_data['y'],
                                                    yr = progress_data['resi'],
                                                    yfit = progress_data['yfit']
                                                )).to_dict('list')

            # model analysis
            if len(list(experiment_df)) > 2:
                model_dict = ck.kinetic_model(experiment_db)
                model_result = model_dict.model(subtract, transform, threshold, bottom, top, slope, scalex, offset)
                model_data_source.data = pd.DataFrame(data = dict(
                                                            xt = model_result['xt'],
                                                            yt = model_result['yt'],
                                                            et = model_result['et'],
                                                            n = model_result['n'],
                                                            ct = model_result['ct']
                                                           )).to_dict('list')
                model_plot_source.data = pd.DataFrame(data = dict(
                                                            xp = model_result['xp'],
                                                            yp = model_result['yp'],
                                                            u = model_result['u'],
                                                            l = model_result['l'],
                                                            ep = model_result['ep'],
                                                            cp = model_result['cp']
                                                           )).to_dict('list')
                model_fit_source.data = pd.DataFrame(data = dict(
                                                            x = model_result['xfit'],
                                                            y = model_result['yfit']
                                                            )).to_dict('list')

                if experiment_db['model'] == 'Michaelis-Menten' or experiment_db['model'] == 'kcat/Km':
                    mm_source.data = pd.DataFrame(data = dict(
                                                           label = ['Fit Value', 'Std. Error'],
                                                           Km = model_result['Km'],
                                                           Vmax = model_result['Vmax'],
                                                           ksp = model_result['ksp'],
                                                           kcat = model_result['kcat']
                                                       ), index=['value', 'error']).to_dict('list')
                    ic_source.data = pd.DataFrame(data = dict(
                                                        label = [], Bottom = [], Top = [],
                                                        Slope = [], p50 = []
                                                    )).to_dict('list')
                    model.xaxis.axis_label = 'Concentration'

                elif experiment_db['model'] == 'pEC50/pIC50':
                    ic_source.data = pd.DataFrame(data = dict(
                                                        label = ['Fit Value', 'Std. Error'],
                                                        Bottom = model_result['Bottom'],
                                                        Top = model_result['Top'],
                                                        Slope = model_result['Slope'],
                                                        p50 = model_result['p50']
                                                    ), index = ['value', 'error']).to_dict('list')
                    mm_source.data = pd.DataFrame(data = dict( label = [], Km = [], Vmax = [], ksp = [], kcat = [] )).to_dict('list')
                    model.xaxis.axis_label = 'Log10(Concentration)'

                else:
                    mm_source.data = pd.DataFrame(data = dict( label = [], Km = [], Vmax = [], ksp = [], kcat = [] )).to_dict('list')
                    ic_source.data = pd.DataFrame(data = dict( label = [], Bottom = [], Top = [],
                                                        Slope = [], p50 =[]
                                                       )).to_dict('list')
                    model.xaxis.axis_label = 'Sample #'
#    except Exception as e:
#        error = str(e)
#        error_page(error, 'Error updating plots due to:')

def load_page(): #experiment_df, experiment_db):
    try:
        ########## bokeh plots ##########

        # general plot tools
        plot_tools = 'wheel_zoom, pan, reset, save, hover'
        hover = HoverTool(tooltips=[("(x,y)", "($x, $y)")])

        # progress curve plots
        global raw_source
        raw_source = ColumnDataSource(data=dict(x=[], y=[], yr=[], yfit=[]))
        global raw
        raw = figure(title="Initial Rate Fit", x_axis_label="Time", y_axis_label="Signal",
                     width=350, height=300, tools=plot_tools)
        raw.circle('x', 'y', size=2, source=raw_source, color='gray',
                    selection_color="black", alpha=0.6, nonselection_alpha=0.2, selection_alpha=0.6)
        raw.line('x', 'yfit', source=raw_source, color='red')
        global warning_source
        warning_source = ColumnDataSource(data=dict(x=[0], y=[0], t=['Please enter transform equation! \nMust convert signal to [substrate] \nin Schnell-Mendoza mode (e.g. via \nx/6.22/0.45/0.001 for sample data). \nNote: this transform may need \nto be inverted through multiplying \nby -1 when analyzing experiments \nthat measure increasing product \nconcentration over time)']))
        global warning
        warning = raw.text(x='x', y='y', text='t', text_font_size='12pt', angle=0, source=warning_source)
        warning.visible = False
        global circles_source
        circles_source = ColumnDataSource(data=dict(x=[-.05, -.05, 1.6, 1.6], y=[0, 0.6, 0, 0.6]))
        global circles
        circles = raw.circle(x='x', y='y', alpha=0., source=circles_source)
        circles.visible = False
        global resi
        resi = figure(title="Initial Rate Fit Residuals", x_axis_label="Time", y_axis_label="Residual",
                     width=700, height=200, tools='wheel_zoom,pan,reset')
        resi.yaxis.formatter = BasicTickFormatter(precision=2, use_scientific=True)
        resi.circle('x', 'yr', size=5, source=raw_source, color='grey', alpha=0.6)

        # model plot for titration experiments
        global model_data_source
        model_data_source = ColumnDataSource(data=dict(xt=[], yt=[], n=[], ct=[], et=[]))
        global model_plot_source
        model_plot_source = ColumnDataSource(data=dict(xp=[], yp=[], l=[], u=[], cp=[], ep=[]))
        global model_fit_source
        model_fit_source = ColumnDataSource(data=dict(x=[], y=[]))
        global varea_source
        varea_source = ColumnDataSource(data=dict(x=[], r1=[], r2=[]))
        global model
        model = figure(title='Model Fit', x_axis_label='Concentration', y_axis_label='Rate', width=350,
                      height=300, tools=plot_tools)
        model.circle('xp', 'yp', size=8, source=model_plot_source, color='cp', alpha=0.6)
        model.add_layout(Whisker(source=model_plot_source, base='xp', upper='u', lower='l'))
        model.line('x', 'y', source=model_fit_source, line_width=3, color='black', alpha=0.8)
        model.varea('x', 'r1', 'r2', source=varea_source, color='grey', alpha=0.3)

        ########## bokeh widgets ##########

        # button for selecting progress curve fitting routine
        global fit_button
        fit_button = RadioButtonGroup(labels=['Maximize Slope Magnitude', 'Linear Fit',
                                        'Logarithmic Fit', 'Schnell-Mendoza'], active=0, width=375)
        fit_button.on_change('active', widget_callback)

        # button for selecting progress curve fitting routine
        global scalex_box
        scalex_box = CheckboxButtonGroup(labels=["Transform X-Axis to Log10 Scale"], active=[])
        scalex_box.on_change('active', widget_callback)

        # dropdown menu for selecting titration experiment model
        global model_select
        model_select = Select(title='Choose Model', value='Michaelis-Menten',
                              options=['Michaelis-Menten', 'kcat/Km', 'pEC50/pIC50', 'High-Throughput Screen'], width=350)
        model_select.on_change('value', widget_callback)

        # dropdown menu for selecting blank sample to subtract from remaining titration samples
        global subtract_select
        subtract_select = Select(title='Select Blank Sample for Subtraction', value='',
                                     options=list(experiment_df)[1:]+[''], width=350)
        subtract_select.on_change('value', widget_callback)

        # dropdown menu for selecting titration sample to plot in current view
        global sample_select
        sample_select = Select(title='Y-Axis Sample', value=list(experiment_df)[-1],
                               options=list(experiment_df)[1:], width=350)
        sample_select.on_change('value', sample_callback)

        # text input box for transforming slopes to rates
        global transform_input
        transform_input = TextInput(value='', title="Enter Transform Equation", width=350)
        transform_input.on_change('value', widget_callback)

        # text input box for setting delay time in logarithmic progress curve fitting
        global offset_input
        offset_input = TextInput(value='', title="Enter Time Between Mixing and First Read", width=350)
        offset_input.on_change('value', widget_callback)

        # text input boxes for fixing EC50/IC50 parameters
        global bottom_fix
        bottom_fix = TextInput(value='', title="Fix pEC50/pIC50 Bottom")
        bottom_fix.on_change('value', widget_callback)

        global top_fix
        top_fix = TextInput(value='', title="Fix pEC50/pIC50 Top")
        top_fix.on_change('value', widget_callback)

        global slope_fix
        slope_fix = TextInput(value='', title="Fix pEC50/pIC50 Hill Slope")
        slope_fix.on_change('value', widget_callback)

        # text input boxes for progress curve xrange selection
        global start_time
        start_time = TextInput(value=str(experiment_df[list(experiment_df)[0]].values[0]), title="Enter Start Time")
        global end_time
        end_time = TextInput(value=str(experiment_df[list(experiment_df)[0]].values[-1]), title='Enter End Time')
        start_time.on_change('value', xbox_callback)
        end_time.on_change('value', xbox_callback)

        # range slider to select threshold for hit detection in HTS mode
        global threshold_slider
        threshold_slider = Slider(start=0, end=5, value=2, step=0.1,
                        title='HTS Hit Threshold (Standard Deviation)', width=350)
        threshold_slider.on_change('value', threshold_callback)

        # range slider to update plots according to progress cuve xrange selection
        xmin = experiment_df[experiment_df.columns[0]].values[0]
        xmax = experiment_df[experiment_df.columns[0]].values[-1]
        global range_slider
        range_slider = RangeSlider(start=xmin, end=xmax, value=(xmin, xmax),
                        step=experiment_df[experiment_df.columns[0]].values[1]-xmin,
                        title='Fine Tune X-Axis Range', width=650)
        range_slider.on_change('value', slider_callback)

        # button to upload local data file
        global file_source
        file_source = ColumnDataSource(data=dict(file_contents = [], file_name = []))
        file_source.on_change('data', file_callback)
        try:
            output_filename = file_source.data['file_name']+'-out.csv'
        except:
            output_filename = 'output.csv'
        global upload_button
        upload_button = Button(label="Upload Local File", button_type="success", width=350)
        upload_button.js_on_event(events.ButtonClick, CustomJS(args=dict(file_source=file_source),
                                   code=open(join(dirname(__file__), "upload.js")).read()))

        # table containing rate fits and errors
        template="""
        <div style="background:<%=ct%>"; color="white";>
        <%= value %></div>
        """
        formatter = HTMLTemplateFormatter(template=template)
        columns = [
            TableColumn(field='n', title='Sample'),
            TableColumn(field='yt', title='Slope (Initial Rate)', formatter=formatter),
            TableColumn(field='et', title='Std. Error')
        ]
        global rate_table
        rate_table = DataTable(source=model_data_source, columns=columns, width=350, height=250,
                               selectable=True, editable=True)

        # tables containing model fits and errors
        global mm_source
        mm_source = ColumnDataSource(dict(label = [], Km = [], Vmax = []))
        columns = [
            TableColumn(field='label', title=''),
            TableColumn(field='Vmax', title='Vmax'),
            TableColumn(field='Km', title='Km'),
            TableColumn(field='ksp', title='kcat/Km'),
            TableColumn(field='kcat', title='kcat')
        ]
        global mm_table
        mm_table = DataTable(source=mm_source, columns=columns, width=350, height=75,
                               selectable=True, editable=True)
        global ic_source
        ic_source = ColumnDataSource(dict(label = [], Bottom = [], Top = [], Slope = [], p50 = []))
        columns = [
            TableColumn(field='label', title=''),
            TableColumn(field='Bottom', title='Bottom'),
            TableColumn(field='Top', title='Top'),
            TableColumn(field='Slope', title='Slope'),
            TableColumn(field='p50', title='pEC50/IC50')
        ]
        global ic_table
        ic_table = DataTable(source=ic_source, columns=columns, width=350, height=75,
                               selectable=True, editable=True)

        # button for copying rate data table to clipboard
        global copy_button
        copy_button = Button(label="Copy Table to Clipboard", button_type="primary", width=350)
        copy_button.js_on_event(events.ButtonClick, CustomJS(args=dict(source=model_data_source),
                                   code=open(join(dirname(__file__), "copy.js")).read()))

        # button for downloading rate data table to local csv file
        global download_button
        download_button = Button(label="Download Table to CSV", button_type="primary", width=350)
        download_button.js_on_event(events.ButtonClick, CustomJS(args=dict(source=model_data_source, file_name=output_filename),
                                   code=open(join(dirname(__file__), "download.js")).read()))

        ########## document formatting #########

        desc = Div(text=open(join(dirname(__file__), "description.html")).read(), width=1400)

        advanced = Div(text="""<strong>Advanced Settings for \npEC50/IC50 Analysis</strong>""")

        widgets = column(model_select, sample_select, subtract_select,
                            transform_input, offset_input, advanced, scalex_box, bottom_fix, top_fix, slope_fix)
        table = column(rate_table)
        main_row = row(column(upload_button, widgets),
                        column(fit_button, row(raw, model), resi, row(start_time, end_time), range_slider),
                        column(download_button, copy_button, table, mm_table, ic_table, threshold_slider))

        sizing_mode = 'scale_width'
        l = layout([
            [desc],
            [main_row]
        ], sizing_mode=sizing_mode)

        update()
        curdoc().clear()
        curdoc().add_root(l)
        curdoc().title = "ICEKAT"
    except Exception as e:
        error = str(e)
        error_page(error, 'Error creating page due to:')

def file_callback(attrname, old, new):
        loading = Div(text=open(join(dirname(__file__), "loader.html")).read(), width=1400)
        l = layout([loading], sizing_mode='scale_both')
        curdoc().clear()
        curdoc().add_root(l)
        curdoc().title = "ICEKAT"
        curdoc().add_next_tick_callback(file_callback2)

def file_callback2():
    try:
        # decode data
        raw_contents = file_source.data['file_contents'][0]
        prefix, b64_contents = raw_contents.split(',', 1)
        file_contents = base64.b64decode(b64_contents).decode('utf-8-sig', errors='ignore')
        file_io = StringIO(file_contents)

        # update dataframe
        global experiment_df
        experiment_df = pd.read_csv(file_io)
        experiment_df.columns = [str(i) for i in list(experiment_df)]
        experiment_df = experiment_df.dropna(axis=1, how='all').dropna(axis=0)
        experiment_df = experiment_df.drop_duplicates(experiment_df.columns[0])

        # update database
        global experiment_db
        experiment_db = dict(model = 'Michaelis-Menten')
        xmin = experiment_df[experiment_df.columns[0]].values[0]
        xmax = experiment_df[experiment_df.columns[0]].values[-1]
        for s in experiment_df.columns[1:]:
            #if np.max(experiment_df[s]) > 0.:
            df = experiment_df[[experiment_df.columns[0], s]]
            experiment_db[s] = ck.progress_curve(df, xmin, xmax)
            experiment_db[s+'_fit'] = 0
        load_page()
    except Exception as e:
        error = str(e)
        error_page(error, 'Error loading data file due to:')

def error_page(error, f):
    function = Div(text=f, width=1400)
    error_message = Div(text=error, width=1400)
    reload = Div(text='Please refresh page to try again.', width=1400)
    spacer = Div(text='', width=1400)
    l = layout([function], [error_message], [spacer], [reload], sizing_mode='stretch_both')
    curdoc().clear()
    curdoc().add_root(l)
    curdoc().title = "ICEKAT"

########## import sample data ##########

experiment_file = join(dirname(__file__), 'test.csv')
experiment_df = pd.read_csv(experiment_file)
experiment_df.columns = [str(i) for i in list(experiment_df)]
experiment_db = dict(model = 'Michaelis-Menten')
xmin = experiment_df[experiment_df.columns[0]].values[0]
xmax = experiment_df[experiment_df.columns[0]].values[-1]
for s in experiment_df.columns[1:]:
    df = experiment_df[[experiment_df.columns[0], s]]
    experiment_db[s] = ck.progress_curve(df, xmin, xmax)
    experiment_db[s+'_fit'] = 0

# load page
load_page()
