import numpy as np
from bokeh.layouts import widgetbox, column, row
from bokeh.models import Slider
import bokeh.plotting as plt
from bokeh.models import HoverTool, ColumnDataSource, Select, Slider, WheelZoomTool, Range1d, markers, Arrow, NormalHead
from bokeh.io import output_notebook, show, push_notebook, output_file, curdoc
output_file("test.html")


def radii(gamma, x0, y0, r):
    '''calculates the two intersections of a line from the origin with the circle defining the window'''
    r1 = np.cos(gamma)*x0+np.sin(gamma)*y0 + np.sqrt((np.cos(gamma)*x0+np.sin(gamma)*y0)**2  -(x0**2+y0**2-r**2))
    r2 = np.cos(gamma)*x0+np.sin(gamma)*y0 - np.sqrt((np.cos(gamma)*x0+np.sin(gamma)*y0)**2  -(x0**2+y0**2-r**2))
    return (r1,r2)
def excentricity(beta, e_window, e_sample=[0,-46]):
    '''calculates the excentricity of the circle defining the window for different top flange rotations'''
    x0 = np.cos(beta)*e_window - e_sample[0]
    y0 = np.sin(beta)*e_window - e_sample[1]
    return x0, y0

def phi(r, h):
    '''
    Returns the scatter angle measured from the vertical of r with height above the sample,h.
    '''
    return np.arctan(r/h)

def delta_from_beta(beta, gamma, r, e_window, e_sample, h):
    x0,y0 = excentricity(beta, e_window, e_sample)
    radiis = radii(gamma, x0, y0,r)
    deltas = [90-np.rad2deg(phi(radiius, h)) for radiius in radiis]
    return deltas
    
def circle(beta, e_window, e_sample, r):
    '''calculates the excentricity of the circle defining the window for different top flange rotations'''
    x0,y0 = excentricity(beta, e_window, e_sample)
    phi = np.linspace(0, 2*np.pi, 3601, endpoint=True)
    x = r * np.cos(phi) - x0
    y = r * np.sin(phi) - y0
    return x, y

'''Be window - scattering angles'''
gamma = np.linspace(0, np.pi, 3601, endpoint=True)
h = 36
r = 70/2
e_window = 30
e_sample = [0,-46]
deltamin, deltamax = delta_from_beta(0, gamma, r, e_window, e_sample, h)
source_window = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y1=deltamin, y2 = deltamax))
#source = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y=deltamin))
x,y = circle(0, e_window, e_sample, r)
source_window_circle = ColumnDataSource(data=dict(x=x, y=y))


'''Top flange - scattering angles'''
h_t = 102
r_t = 248.5/2
e_window_t = 0
deltamin_t, deltamax_t = delta_from_beta(0, gamma, r_t, e_window_t, e_sample, h_t)
source_top_flange = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y1=deltamin_t, y2 = deltamax_t))
x,y = circle(0, e_window_t, e_sample, r_t)
source_top_flange_circle = ColumnDataSource(data=dict(x=x, y=y))


g = plt.figure(title='Be Window - Geometry')
imgurl = 'https://drive.switch.ch/index.php/apps/files_sharing/ajax/publicpreview.php?x=3360&y=1126&a=true&file=Chamber_01.jpg&t=ewZ1kPPzxpalNpc&scalingup=0'


g.image_url(url=[imgurl], x=-152, y=-200, w=305, h=305, anchor="bottom_left")
g.line('x', 'y', source=source_window_circle, line_width=5, line_color = 'red', legend_label="Be window")
g.line('x', 'y', source=source_top_flange_circle, line_width=5, line_color = 'red')
#g.image_url(url=['https://docs.bokeh.org/en/latest/_static/images/logo.png'], x=0, y=0, w=100, h=100, anchor="bottom_left")
g.circle(x=0, y=0, size=15,fill_color="green", legend_label="Sample")



def update_plot(attr, old, new):
    beta = new
    #beta = new
    deltamin, deltamax = delta_from_beta(np.deg2rad(beta), gamma, r, e_window, e_sample, h)
    new_source = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y1=deltamin, y2 = deltamax))
    source_window.data = new_source.data
    x,y = circle(np.deg2rad(beta), e_window, e_sample, r)
    new_source = ColumnDataSource(data=dict(x=x, y=y))
    source_window_circle.data = new_source.data

    #p.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
    #push_notebook()

p = plt.figure(title='Be Window - Scatter Angles')
p.line('x', 'y1', source=source_window, line_width=3, line_alpha=0.6)
p.line('x', 'y2', source=source_window, line_width=3, line_alpha=0.6)
p.line('x', 'y1', source=source_top_flange, line_width=3, line_alpha=0.6, line_color='grey')
p.line('x', 'y2', source=source_top_flange, line_width=3, line_alpha=0.6, line_color='grey')
p.add_tools(HoverTool(tooltips=[("(gamma,delta_min)", "($x, $y)")]))
p.add_tools(HoverTool(tooltips=[("(gamma,delta_max)", "($x, $y)")]))

p.toolbar.active_scroll = p.select_one(WheelZoomTool)

p.xaxis.axis_label = 'gamma, deg'
p.yaxis.axis_label = 'delta, deg'
p.background_fill_color = 'beige'
p.background_fill_alpha = 0.2
slider = Slider(start=0., end=360., step=15., value = 0., title='Position of Be window (0° downstream along beam)')
slider.on_change('value',update_plot)
imgs = row(p,g)
layout = column(imgs, slider)
#plt.show(layout, notebook_handle=True)
#plt.show(p)
p.x_range=Range1d(0, 180)
p.y_range=Range1d(0, 180)
curdoc().title = "THC illustrator"
curdoc().add_root(layout)
