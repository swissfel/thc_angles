import numpy as np
from bokeh.layouts import widgetbox, column, row
from bokeh.models import Slider
import bokeh.plotting as plt
from bokeh.models import HoverTool, ColumnDataSource, Select, Slider, WheelZoomTool, Range1d, markers, Arrow, NormalHead
from bokeh.io import output_notebook, show, push_notebook, output_file, curdoc
#output_file("test.html")
from scipy.spatial.transform import Rotation
from scipy.linalg import det

### functions for rectangle intersections
def rotate_detector_2d(l=80, b=40, e_det=10, beta=0):
    """rotates the detector around the center of the flange by beta, where the long detector axis is along x for beta = 0"""
    e1 = np.array([b/2,  l/2+e_det, 0])
    e2 = np.array([-b/2, l/2+e_det, 0])
    e3 = np.array([-b/2, -l/2+e_det, 0])
    e4 = np.array([b/2,  -l/2+e_det, 0])
    r_beta = Rotation.from_euler('z', [beta], degrees=False)
    e1,e2,e3,e4 = r_beta.apply([e1,e2,e3,e4])
    return np.array([e1[:2],e2[:2],e3[:2],e4[:2]])


def radii_rectangle(l=80, b=40, e_det=10, e_sample=np.array([0,46]), gam=np.array([0]), beta=0):
    """calculates the distance between sample and the intersection points of the detector projected onto the x-z pllane for given beta and gamma angles"""
    # x-axis defined in opposite direction
    gam = -gam
    s = e_sample
    s1 = np.zeros((len(gam),2,2)) + s

    s1[:,0,0] = s1[:,0,0] -250*np.array(np.cos(gam))
    s1[:,0,1] = s1[:,0,1] +250*np.array(np.sin(gam))
    s1[:,1,0] = s1[:,1,0] +250*np.array(np.cos(gam))
    s1[:,1,1] = s1[:,1,1] -250*np.array(np.sin(gam))
    
    e1,e2,e3,e4 = rotate_detector_2d(l=l, b=b, e_det=e_det, beta=beta)
    es = np.array([[e1,e2], [e1,e4], [e2,e3], [e3,e4]])
    M = np.array( [np.array([np.divide( det([p1-p3,p3-p4]),det([p1-p2,p3-p4])) for p3, p4 in zip(s1[:,0],s1[:,1]) ]) for p1,p2 in es])
    M[~((M<1.)&(M>0.))] = np.nan
    
    p1s = es[:,0]
    p2s = es[:,1]
    intersects = np.array([m[:,None]*(p2s-p1s)+p1s for m in M.T])
    radii_2d = np.linalg.norm(intersects-s, axis=2)*(intersects-s)[:,:,1]/-abs((intersects-s)[:,:,1])
    radii_2d = np.array([np.nanmax(radii_2d, axis=1), np.nanmin(radii_2d, axis=1)])    
    
    return radii_2d

### functions for circle intersections


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

'''In-vacuum detector - scattering angles 1 position'''
gamma_invac = np.linspace(0, np.pi, 361, endpoint=True)
h_invac = 42
l_invac = 80
b_invac = 40
e_invac = 36
e_sample_invac = np.array([0,46])
beta=np.deg2rad(0)
radii_2d = radii_rectangle(l=l_invac, b=b_invac, e_det=e_invac, e_sample=e_sample_invac, gam=gamma_invac, beta=beta )
deltamin, deltamax = [90-np.rad2deg(phi(radii_2d[0], h_invac)),90-np.rad2deg(phi(radii_2d[1], h_invac))]
source_invac = ColumnDataSource(data=dict(x=np.rad2deg(gamma_invac), y1=deltamin, y2 = deltamax))
rect_edges_invac = rotate_detector_2d(l=l_invac, b=b_invac, e_det=e_invac, beta=beta) - e_sample_invac
rect_edges_invac = np.vstack([rect_edges_invac,rect_edges_invac[0]])
source_invac_rect = ColumnDataSource(data=dict(x=rect_edges_invac.T[0], y=rect_edges_invac.T[1]))

'''In-vacuum detector - scattering angles cumulative all position'''
betas = np.linspace(0,2*np.pi, 120)
gammas = np.linspace(0.001,np.pi+0.001, 120)
radii_2d = np.array([radii_rectangle(l=l_invac, b=b_invac, e_det=e_invac, e_sample=e_sample_invac, gam=gammas, beta=beta ) for beta in betas])
rs = np.array([np.nanmax(radii_2d, axis=(0,1)), np.nanmin(radii_2d, axis=(0,1))])
deltamin, deltamax = [90-np.rad2deg(phi(rs[0], h_invac)),90-np.rad2deg(phi(rs[1], h_invac))]
source_invac_full_betas = ColumnDataSource(data=dict(x=np.rad2deg(gammas), y1=deltamin, y2 = deltamax))



'''Beryllium window - scattering angles 1 position'''
gamma = np.linspace(0, np.pi, 3601, endpoint=True)
h = 36
r = 70/2
e_window = 30
e_sample = [0,-46]
beta=0
deltamin, deltamax = delta_from_beta(np.deg2rad(beta-90), gamma, r, e_window, e_sample, h)
source_window = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y1=deltamin, y2 = deltamax))
#source = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y=deltamin))
x,y = circle(np.deg2rad(beta-90), e_window, e_sample, r)
source_window_circle = ColumnDataSource(data=dict(x=x, y=y))

'''Be window - scattering angles cumulative all position'''
betas = np.linspace(0,2*np.pi, 360)
delta_from_beta(np.deg2rad(0), gamma, r, e_window, e_sample, h)
ds = np.array([delta_from_beta(beta, gamma, r, e_window, e_sample, h) for beta in betas])
dmins =  np.array([np.nanmin(dmin) for dmin in ds[:,0,:].T])
dmaxs = np.array([np.nanmax(dmax) for dmax in ds[:,1,:].T])
source_window_full_betas = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y1=dmins, y2 = dmaxs))

'''Top flange - scattering angles'''
h_t = 102
r_t = 248.5/2
e_window_t = 0
deltamin_t, deltamax_t = delta_from_beta(0, gamma, r_t, e_window_t, e_sample, h_t)
source_top_flange = ColumnDataSource(data=dict(x=np.rad2deg(gamma), y1=deltamin_t, y2 = deltamax_t))
x,y = circle(0, e_window_t, e_sample, r_t)
source_top_flange_circle = ColumnDataSource(data=dict(x=x, y=y))


g = plt.figure(title='Geometry')
imgurl = 'https://drive.switch.ch/index.php/apps/files_sharing/ajax/publicpreview.php?x=3360&y=1040&a=true&file=Be_window9.JPG&t=WgAa3cg9sgVJFZK&scalingup=0'

g.image_url(url=[imgurl], x=-165.5, y=-202.5, w=330, h=350, anchor="bottom_left")
g.line('x', 'y', source=source_window_circle, line_width=5, line_color = 'red', legend_label="Be window")
g.line('x', 'y', source=source_top_flange_circle, line_width=5, line_color = 'red')
g.line('x', 'y', source=source_invac_rect, line_width=5, line_color = 'green', legend_label="In-vacuum detector")
g.circle(x=0, y=0, size=10,fill_color="orange", line_color = 'black', legend_label="Sample")
g.y_range=Range1d(-205, 145)
g.x_range=Range1d(-165, 165)
g.title.align='center'



def update_plot(attr, old, new):
    beta = new-90
    deltamin, deltamax = delta_from_beta(np.deg2rad(beta), gamma, r, e_window, e_sample, h)
    source_window.data = dict(x=np.rad2deg(gamma), y1=deltamin, y2 = deltamax)
    
    radii_2d = radii_rectangle(l=l_invac, b=b_invac, e_det=e_invac, e_sample=e_sample_invac, gam=gamma_invac, beta=np.deg2rad(new) )
    deltamin, deltamax = [90-np.rad2deg(phi(radii_2d[0], h_invac)),90-np.rad2deg(phi(radii_2d[1], h_invac))]
    source_invac.data = dict(x=np.rad2deg(gamma_invac), y1=deltamin, y2 = deltamax)
    rect_edges_invac = rotate_detector_2d(l=l_invac, b=b_invac, e_det=e_invac, beta=np.deg2rad(new)) - e_sample_invac
    rect_edges_invac = np.vstack([rect_edges_invac,rect_edges_invac[0]])
    source_invac_rect.data = dict(x=rect_edges_invac.T[0], y=rect_edges_invac.T[1])
    
    x,y = circle(np.deg2rad(beta), e_window, e_sample, r)
    source_window_circle.data = dict(x=x, y=y)

    #p.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
    #push_notebook()

p = plt.figure(title='Diffraction Angles')
p.line('x', 'y1', source=source_window, line_width=3, line_alpha=0.6, legend_label="Be window current position")
p.line('x', 'y2', source=source_window, line_width=3, line_alpha=0.6)
p.line('x', 'y1', source=source_window_full_betas, line_width=3, line_alpha=0.6, line_dash='dashed')
p.line('x', 'y1', source=source_invac, line_width=3, line_alpha=0.6, line_color='green', legend_label="In-vacuum detector current position")
p.line('x', 'y2', source=source_invac, line_width=3, line_alpha=0.6, line_color='green')
p.line('x', 'y1', source=source_invac_full_betas, line_width=3, line_alpha=0.6, line_color='green', line_dash='dashed')
p.line('x', 'y2', source=source_invac_full_betas, line_width=3, line_alpha=0.6, line_color='green', line_dash='dashed')
p.line('x', 'y2', source=source_window_full_betas, line_width=3, line_alpha=0.6, line_dash='dashed')
p.line('x', 'y1', source=source_top_flange, line_width=3, line_alpha=0.6, line_color='grey', legend_label="Top of CF250 flange")
p.line('x', 'y2', source=source_top_flange, line_width=3, line_alpha=0.6, line_color='grey')

p.add_tools(HoverTool(tooltips=[("(gamma,delta)", "($x, $y)")]))

p.toolbar.active_scroll = p.select_one(WheelZoomTool)

p.xaxis.axis_label = 'gamma, deg'
p.yaxis.axis_label = 'delta, deg'
#p.background_fill_color = 'beige'
#p.background_fill_alpha = 0.2
p.legend.location='top_center'
p.title.align='center'
slider = Slider(start=0., end=360., step=15., value = 0., title='Position of Be window / in-vacuum detector (90° downstream along beam)')
slider.on_change('value',update_plot)
imgs = row(p,g)
layout = column(imgs, slider)
#plt.show(layout, notebook_handle=True)
#plt.show(p)
p.x_range=Range1d(0, 180)
p.y_range=Range1d(0, 180)
curdoc().title = "THC illustrator"
curdoc().add_root(layout)
