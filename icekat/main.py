import continuous_kinetics as ck
from os.path import dirname, join
from io import StringIO
import base64
import pandas as pd
import numpy as np
from bokeh.io import curdoc
from bokeh.layouts import row, column, layout
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, Div, BasicTickFormatter, Whisker, Column
from bokeh.models.widgets import RadioButtonGroup, Select, TextInput, Button, DataTable, TableColumn, RangeSlider, Slider, HTMLTemplateFormatter, CheckboxButtonGroup
from bokeh.plotting import figure

########## bokeh methods ##########
def widget_callback(attrname, old, new):

    start_e = float(start_time.value)
    stop_e = float(end_time.value)
    range_slider.value = (start_e, stop_e)

    raw_title = raw.title.text
    resi_title = resi.title.text
    model_title = model.title.text
    raw.title.text = "Calculating..."
    resi.title.text = "Calculating..."
    model.title.text = "Calculating..."

    fit_routine = fit_button.active
    for s in experiment_df.columns[1:]:
        df = experiment_df[[experiment_df.columns[0], s]]
        xmin, xmax = min(df[df.columns[0]]), max(df[df.columns[0]])
        experiment_db[s] = ck.progress_curve(df, xmin, xmax)
        if fit_routine == 0:
            progress_data = experiment_db[s].spline()
            experiment_db[s].spline = progress_data
            experiment_db[s+'_fit'] = 0
        elif fit_routine == 1:
            progress_data = experiment_db[s].linear()
            experiment_db[s].linear = progress_data
            experiment_db[s+'_fit'] = 1
        #elif fit_routine == 2:
        #    try:
        #        offset = float(offset)
        #    except:
        #        offset = 0.0
        #    progress_data = experiment_db[s].logarithmic(offset)
        #    experiment_db[s].logarithmic = progress_data
        #    experiment_db[s+'_fit'] = 2

    raw.title.text = raw_title
    resi.title.text = resi_title
    model.title.text= model_title

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

    if start >= stop -  5*step and stop <= end - 5*step:
        range_slider.value = (start, start + 5*step)

    elif start >= stop - 5*step and stop > end - 5*step:
        range_slider.value = (stop - 5*step, stop)

    if start_e != start or stop_e != stop:
        start_time.value = str(start)
        end_time.value = str(stop)

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

    # update database
    experiment_db[sample+'_fit'] = fit_routine
    experiment_db['model'] = model_eq
    pdf = experiment_df[[experiment_df.columns[0], sample]]
    experiment_db[sample] = ck.progress_curve(pdf, start, end)

    if fit_routine == 3:
        raw.title.text = "Schnell-Mendoza Fit"
        resi.title.text = "Schnell-Mendoza Fit Residuals"
        model_select.value = 'Michaelis-Menten'
        scalex_box.active = []
        raw_data, model_result, fit_data, varea_data = ck.sm_fit(experiment_db).fit(sample, transform, subtract)
        raw_source.data = pd.DataFrame(data = dict(
                                             x = raw_data['x'], y = raw_data['y'],
                                             yr = raw_data['resi'], yfit = raw_data['yfit'],
                                            )).to_dict('list')
        model_data_source.data = pd.DataFrame(data = dict(
                                                    xt = [], yt = [], et = [],
                                                    n = [], ct = []
                                                    )).to_dict('list')
        model_plot_source.data = pd.DataFrame(data = dict(
                                                    xp = [], yp = [], u = [],
                                                    l = [], ep = [], cp = []
                                                    )).to_dict('list')
        model_fit_source.data = pd.DataFrame(data = dict(
                                                   x = model_result['xfit'],
                                                   y = model_result['yfit']
                                                   )).to_dict('list')
        varea_source.data = pd.DataFrame(data = dict(
                                               x = varea_data['x'], r1 = varea_data['r1'],
                                               r2 = varea_data['r2']
                                              )).to_dict('list')
        mm_source.data = pd.DataFrame(data = dict(
                                             label = ['Fit Value', 'Std. Error'],
                                             Km = fit_data['Km'],
                                             Vmax = fit_data['Vmax']
                                           ), index=['value', 'error']).to_dict('list')
        ic_source.data = pd.DataFrame(data = dict(
                                             label = [], Bottom = [], Top = [],
                                             Slope = [], p50 = []
                                           )).to_dict('list')
        model.xaxis.axis_label = 'Concentration'
        raw.yaxis.axis_label = 'Concentration'

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

            if experiment_db['model'] == 'Michaelis-Menten':
                mm_source.data = pd.DataFrame(data = dict(
                                                       label = ['Fit Value', 'Std. Error'],
                                                       Km = model_result['Km'],
                                                       Vmax = model_result['Vmax']
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
                mm_source.data = pd.DataFrame(data = dict( label = [], Km = [], Vmax = [] )).to_dict('list')
                model.xaxis.axis_label = 'Log10(Concentration)'

            else:
                mm_source.data = pd.DataFrame(data = dict( label = [], Km = [], Vmax = [] )).to_dict('list')
                ic_source.data = pd.DataFrame(data = dict( label = [], Bottom = [], Top = [],
                                                    Slope = [], p50 =[]
                                                   )).to_dict('list')
                model.xaxis.axis_label = 'Sample #'

def load_page(experiment_df, experiment_db):

    ########## bokeh plots ##########

    # general plot tools
    plot_tools = 'wheel_zoom, pan, reset, save, hover'
    hover = HoverTool(tooltips=[("(x,y)", "($x, $y)")])

    # progress curve plots
    global raw_source
    raw_source = ColumnDataSource(data=dict(x=[], y=[], yr=[], yfit=[]))
    global raw
    raw = figure(title="Initial Rate Fit", x_axis_label="Time", y_axis_label="Signal",
                 plot_width=350, plot_height=300, tools=plot_tools)
    raw.circle('x', 'y', size=2, source=raw_source, color='gray',
                selection_color="black", alpha=0.6, nonselection_alpha=0.2, selection_alpha=0.6)
    raw.line('x', 'yfit', source=raw_source, color='red')
    global resi
    resi = figure(title="Initial Rate Fit Residuals", x_axis_label="Time", y_axis_label="Residual",
                 plot_width=700, plot_height=200, tools='wheel_zoom,pan,reset')
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
    model = figure(title='Model Fit', x_axis_label='Concentration', y_axis_label='Rate', plot_width=350,
                  plot_height=300, tools=plot_tools)
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
    scalex_box = CheckboxButtonGroup(labels=["transform x-axis to Log10 scale"], active=[])
    scalex_box.on_change('active', widget_callback)

    # dropdown menu for selecting titration experiment model
    global model_select
    model_select = Select(title='Choose Model', value='Michaelis-Menten',
                          options=['Michaelis-Menten', 'pEC50/pIC50', 'High-Throughput Screen'], width=350)
    model_select.on_change('value', widget_callback)

    # dropdown menu for selecting blank sample to subtract from remaining titration samples
    global subtract_select
    subtract_select = Select(title='Select Blank Sample for Subtraction', value='',
                                 options=list(experiment_df)[1:]+[''], width=350)
    subtract_select.on_change('value', widget_callback)

    # dropdown menu for selecting titration sample to plot in current view
    global sample_select
    sample_select = Select(title='Y Axis Sample', value=list(experiment_df)[-1],
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

    # text input boxes for fixing EC/IC50 parameters
    global bottom_fix
    bottom_fix = TextInput(value='', title="Fix pIC50/pEC50 Bottom")
    bottom_fix.on_change('value', widget_callback)

    global top_fix
    top_fix = TextInput(value='', title="Fix pIC50/pEC50 Top")
    top_fix.on_change('value', widget_callback)

    global slope_fix
    slope_fix = TextInput(value='', title="Fix pIC50/pEC50 Hill Slope")
    slope_fix.on_change('value', widget_callback)

    # text input boxes for progress curve xrange selection
    global start_time
    start_time = TextInput(value=str(experiment_df[list(experiment_df)[0]].values[0]), title="Start Time")
    global end_time
    end_time = TextInput(value=str(experiment_df[list(experiment_df)[0]].values[-1]), title='End Time')
    start_time.on_change('value', widget_callback)
    end_time.on_change('value', widget_callback)

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
                    title='X-Axis Range', width=650)
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
    callback = CustomJS(args=dict(file_source=file_source), code=open(join(dirname(__file__), "upload.js")).read())
    upload_button.js_on_event('tap', callback)

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
        TableColumn(field='Km', title='Km')
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
        TableColumn(field='p50', title='pEC/IC50')
    ]
    global ic_table
    ic_table = DataTable(source=ic_source, columns=columns, width=350, height=75,
                           selectable=True, editable=True)

    # button for copying rate data table to clipboard
    global copy_button
    copy_button = Button(label="Copy Table to Clipboard", button_type="primary", width=350)
    callback = CustomJS(args=dict(source=model_data_source),
                               code=open(join(dirname(__file__), "copy.js")).read())
    copy_button.js_on_event('tap', callback)

    # button for downloading rate data table to local csv file
    global download_button
    download_button = Button(label="Download Table to CSV", button_type="primary", width=350)
    callback = CustomJS(args=dict(source=model_data_source, file_name=output_filename),
                               code=open(join(dirname(__file__), "download.js")).read())
    download_button.js_on_event('click', callback)

    ########## document formatting #########

    desc = Div(text=open(join(dirname(__file__), "description.html")).read(), width=1400)

    advanced = Div(text="""<strong>Advanced Settings for \npEC/IC50 Analysis</strong>""")

    widgets = Column(model_select, sample_select, subtract_select,
                        transform_input, offset_input, advanced, scalex_box, bottom_fix, top_fix, slope_fix)
    table = Column(rate_table)
    main_row = row(column(upload_button, widgets),
                    column(fit_button, row(raw, model), resi, range_slider, row(start_time, end_time)),
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

def file_callback(attrname, old, new):

    # decode data
    raw_contents = file_source.data['file_contents'][0]
    prefix, b64_contents = raw_contents.split(',', 1)
    file_contents = base64.b64decode(b64_contents).decode('utf-8-sig', errors='ignore')
    file_io = StringIO(file_contents)

    # update dataframe
    global experiment_df
    experiment_df = pd.read_csv(file_io)
    experiment_df.columns = [str(i) for i in list(experiment_df)]

    # update database
    global experiment_db
    experiment_db = dict(model = 'Michaelis-Menten')
    xmin = experiment_df[experiment_df.columns[0]].values[0]
    xmax = experiment_df[experiment_df.columns[0]].values[-1]
    for s in experiment_df.columns[1:]:
        if np.max(experiment_df[s]) > 0.:
            df = experiment_df[[experiment_df.columns[0], s]].dropna(axis=0)
            experiment_db[s] = ck.progress_curve(df, xmin, xmax)
            experiment_db[s+'_fit'] = 0

    # reload page
    load_page(experiment_df, experiment_db)

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
load_page(experiment_df, experiment_db)
