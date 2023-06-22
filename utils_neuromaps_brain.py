import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from neuromaps.datasets import fetch_fslr
from surfplot import Plot
from surfplot.datasets import load_example_data
from neuromaps.datasets import fetch_fslr
from neuromaps.transforms import mni152_to_fslr
from brainspace.datasets import load_parcellation
import matplotlib.ticker as ticker

import os
import numpy as np
import svgutils.transform as sg
import copy


def plot_adjustments(ax=False, ticks_width=1.5, axis_width=2.5):
    if ax == False:
        ax = plt.gca()
    plt.rcParams['font.family'] = "PT Serif Caption"
    plt.rcParams['xtick.major.width'] = ticks_width
    plt.rcParams['ytick.major.width'] = ticks_width
    plt.rcParams['axes.linewidth'] = ticks_width
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_linewidth(axis_width)
    ax.spines["right"].set_linewidth(axis_width)
    ax.spines["left"].set_linewidth(axis_width)
    ax.spines["bottom"].set_linewidth(axis_width)
    ax.tick_params(direction='out', length=4, width=axis_width + 0.5, colors='k',
                   grid_color='k', grid_alpha=0.5, axis='both')
    ax.tick_params(direction='out', which='minor', length=6, width=axis_width + 0.5, colors='k',
                   grid_color='k', grid_alpha=0.5, axis='both')
    ax.tick_params(direction='out', which='major', length=6, colors='k',
                   grid_color='k', grid_alpha=0.5, axis='both')
    ax.tick_params(axis='both', which='major', labelsize=16)


# Truncating a colormap (between two intervals)
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    '''
    https://stackoverflow.com/a/18926541
    '''
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def normal_view(current_nodestrength, edges=False, surftype='inflated',
                xlabel=r'$<s_i> $', q_thresh=0.0, cmap='RdBu_r',
                brightness=0.7, exp_form=True, parcellation=100,
                vmin=None, vmax=None,graymap_rev=False,alpha_graymap=1.,parcellation_name='schaefer',center_cbar=False):
    '''
    Compute 2 different brain views to obtain a simple representation of the brain map
    by considering 'lateral' and 'medial'views

    Parameters
    ----------
    current_nodestrength: vectors of real 
            Vector containing the real values to be plotted
            on the brain surface
    edges: Bool, True or False, optional
            Activate/disactivate the countour of the parcellation
            in gray (default: False) 
    surftype: string, either 'inflated','veryinflated' or 'midthickness'
            Type of surface selected for rendering (default: 'inflated').
    xlabel: string, optional
            Label used for the legend (default: 's_i')
    q_thresh: real between 0 and 1, optional
            Percentile used to threshold the weights, the corresponding
            regions lower than q_thresh will be displayed in gray
    cmap: colormap, optional
            Color map used for the cortical surface (default: 'custom',
            which is obtained taking the `upper` part of 'jet')
    brightness: real between 0 and 1, optional
            Level of brightness of the surface (default: 0.7). Decrease
            its value to obtain darker surfaces
    exp_form: Bool, True or False, optional
            Enable or disable the exponential notation for the ticks
            of the legend (default: True)
    '''
    if cmap == 'custom':
        cmap_base = 'jet'
        vmincolor, vmaxcolor = 0.5, 1
        cmap = truncate_colormap(cmap_base, vmincolor, vmaxcolor)
    elif cmap == 'RdBu_r':
        cmap='RdBu_r'
    
    mygraymap=mpl.cm.gray
    if graymap_rev==True:
        mygraymap=mpl.cm.gray_r
        
    if alpha_graymap!=1:
        tempcmap = mygraymap
        # Get the colormap colors
        my_cmap = tempcmap(np.arange(tempcmap.N))

        # Set alpha
        my_cmap[:,-1] = [alpha_graymap]*len(my_cmap[:,-1])

        # Create new colormap
        my_cmap = ListedColormap(my_cmap)
        mygraymap=my_cmap
    
        
    # Choose colormap
#     cmap = mpl.cm.jet

#     # Get the colormap colors
#     my_cmap = cmap(np.arange(cmap.N))

#     # Set alpha
#     my_cmap[:,-1] = np.linspace(0.5, 1, cmap.N)
# #     my_cmap[:,-1] = [0.5]*len(my_cmap[:,-1])

#     # Create new colormap
#     my_cmap = ListedColormap(my_cmap)
#     cmap=my_cmap

    surfaces = fetch_fslr(density='32k')
    lh, rh = surfaces[surftype]
    label = xlabel
    p = Plot(lh, rh, zoom=1.5, flip=False,
             brightness=brightness, size=(800, 600))
    lh_sulc, rh_sulc = surfaces['sulc']
    lh_parc, rh_parc = load_parcellation(parcellation_name, scale=parcellation)
    # Modify the value of the parcellation
    lh_parc_mod = np.zeros(len(lh_parc))
    rh_parc_mod = np.zeros(len(rh_parc))
    for idx, l in enumerate(lh_parc):
        if l != 0:  # only for values different zero (zero is subcortical)
            lh_parc_mod[idx] = current_nodestrength[l - 1]
    for idx, l in enumerate(rh_parc):
        if l != 0:  # only for values different zero (zero is subcortical)
            rh_parc_mod[idx] = current_nodestrength[l - 1]
        # Find regions with a value greater that q_thresh quartile
    thresh = np.quantile(current_nodestrength, q=q_thresh)
    lh_regions = np.where(lh_parc_mod >= thresh, lh_parc_mod, 0)
    rh_regions = np.where(rh_parc_mod >= thresh, rh_parc_mod, 0)
    
    
    

    # Do all the lateral
    if vmin == None and vmax == None:
        if q_thresh==0.0:
            p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=True, cmap=cmap,
                        color_range=(min(current_nodestrength),max(current_nodestrength)))
        if q_thresh!=0.0:
            p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=True, cmap=cmap,
                        color_range=(thresh, max(current_nodestrength)))
    else:
        if q_thresh==0.0:
            p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=True, cmap=cmap,
                        color_range=(vmin,vmax))
        if q_thresh!=0.0:
            p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=True, cmap=cmap,
                        color_range=(thresh, vmax))
        
                        
    if edges == True:
        lh_regions_gray=np.random.normal(size=len(lh_regions))+np.mean(lh_regions)
        rh_regions_gray=np.random.normal(size=len(lh_regions))
        p.add_layer({'left': lh_regions, 'right': rh_regions}, cmap=mygraymap,
                    as_outline=True, cbar=False)
    if center_cbar == True:
        pad = -0.51
    else:
        pad = 0.08
        
    
    kws = {'location': 'bottom', 'pad':pad,'label_direction': 0, 'decimals': 1,
           'fontsize': 16, 'n_ticks': 2, 'shrink': .25, 'aspect': 12,
           'draw_border': True}
    fig = p.build(scale=(2, 2), cbar_kws=kws)
    plot_adjustments()
    fig.axes[1].set_xlabel(label, labelpad=-12,
                           weight='normal', fontstyle='normal', fontsize=20)
    
    
    # Adjusting the format of the colorbar range
    if exp_form == True:
        xmin_cbar, xmax_cbar = fig.axes[1].get_xlim()
        scale_pow = int(np.log10(xmax_cbar))

        def my_formatter_fun(x, p):
            return "%.2f" % (x / (10 ** scale_pow))
        fig.axes[1].get_xaxis().set_major_formatter(
            ticker.FuncFormatter(my_formatter_fun))
        fig.axes[1].set_xlabel(label, fontsize=20)
        if center_cbar == True:
            annotate_pad=-pad+0.19
        else:
            annotate_pad=pad
        plt.annotate('x $10^{{{0:d}}}$'.format(scale_pow), xy=(0.6, -0.16+annotate_pad),
                     fontsize=14, xycoords='axes fraction', color='k', alpha=0.6)
    # plt.close()
    return(fig)


def full_view(current_nodestrength, edges=False, surftype='inflated', xlabel=r'$<s_i> $', q_thresh=0.25, cmap='custom',
    brightness=0.7, exp_form=True, parcellation=100, vmin_plot=None, vmax_plot=None,graymap_rev=False,alpha_graymap=1.,parcellation_name='schaefer'):
    '''
    Compute 6 different brain views to obtain a full representation of the brain map
    by considering 'lateral','medial','dorsal','ventral','anterior','posterior' views

    Parameters
    ----------
    current_nodestrength: vectors of real 
            Vector containing the real values to be plotted
            on the brain surface
    edges: Bool, True or False, optional
            Activate/disactivate the countour of the parcellation
            in gray (default: False) 
    surftype: string, either 'inflated','veryinflated' or 'midthickness'
            Type of surface selected for rendering (default: 'inflated').
    xlabel: string, optional
            Label used for the legend (default: 's_i')
    q_thresh: real between 0 and 1, optional
            Percentile used to threshold the weights, the corresponding
            regions lower than q_thresh will be displayed in gray
    cmap: colormap, optional
            Color map used for the cortical surface (default: 'custom',
            which is obtained taking the `upper` part of 'jet')
    brightness: real between 0 and 1, optional
            Level of brightness of the surface (default: 0.7). Decrease
            its value to obtain darker surfaces
    exp_form: Bool, True or False, optional
            Enable or disable the exponential notation for the ticks
            of the legend (default: True)
    '''
    if cmap == 'custom':
        cmap_base = 'jet'
        vmin, vmax = 0.5, 1
        cmap = truncate_colormap(cmap_base, vmin, vmax)

    mygraymap=mpl.cm.gray
    if graymap_rev==True:
        mygraymap=mpl.cm.gray_r
        
    if alpha_graymap!=1:
        tempcmap = mygraymap
        # Get the colormap colors
        my_cmap = tempcmap(np.arange(tempcmap.N))

        # Set alpha
        my_cmap[:,-1] = [alpha_graymap]*len(my_cmap[:,-1])

        # Create new colormap
        my_cmap = ListedColormap(my_cmap)
        mygraymap=my_cmap

    surfaces = fetch_fslr(density='32k')
    lh, rh = surfaces[surftype]
    label = xlabel
    p = Plot(lh, rh, zoom=1.22, layout='grid', views=[
             'lateral'], flip=False, brightness=brightness, size=(2200, 600))
    lh_sulc, rh_sulc = surfaces['sulc']
    lh_parc, rh_parc = load_parcellation(parcellation_name, scale=parcellation)
    # Modify the value of the parcellation
    lh_parc_mod = np.zeros(len(lh_parc))
    rh_parc_mod = np.zeros(len(rh_parc))
    for idx, l in enumerate(lh_parc):
        if l != 0:  # only for values different zero (zero is subcortical)
            lh_parc_mod[idx] = current_nodestrength[l-1]
    for idx, l in enumerate(rh_parc):
        if l != 0:  # only for values different zero (zero is subcortical)
            rh_parc_mod[idx] = current_nodestrength[l-1]
        # Find regions with a value greater that q_thresh quartile
    thresh = np.quantile(current_nodestrength, q=q_thresh)
    lh_regions = np.where(lh_parc_mod >= thresh, lh_parc_mod, 0)
    rh_regions = np.where(rh_parc_mod >= thresh, rh_parc_mod, 0)
    
    
    upper_thresh = max(current_nodestrength)
    if vmin_plot != None and vmax_plot != None:
        thresh = vmin_plot
        upper_thresh = vmax_plot

    # Do all the lateral
#     if vmin==None and vmax==None:
    p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=False, cmap=cmap,
                color_range=(thresh, upper_thresh))
#     else:
#         p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=True,cmap=cmap,
#                 color_range=(vmin,vmax))
    if edges == True:
        p.add_layer({'left': lh_regions, 'right': rh_regions}, cmap=mygraymap,
                    as_outline=True, cbar=False)
    fig = p.build(scale=(2, 2))
    fig.savefig('all_lateral.svg', dpi=150, transparent=True)
    plt.close()

    # Do the the medial view
    p = Plot(lh, rh, zoom=1.22, layout='row', views=[
             'medial'], flip=False, brightness=brightness, size=(2300, 600))
    p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=True, cmap=cmap,
                color_range=(thresh, upper_thresh))
    if edges == True:
        p.add_layer({'left': lh_regions, 'right': rh_regions}, cmap=mygraymap,
                    as_outline=True, cbar=False)
    kws = {'location': 'bottom', 'label_direction': 0, 'decimals': 2,
           'fontsize': 26, 'n_ticks': 2, 'shrink': .25, 'aspect': 12,
           'draw_border': True, 'pad': 0.21}
    fig = p.build(scale=(2, 2), cbar_kws=kws)
    plot_adjustments()
    fig.axes[1].set_xlabel(label, labelpad=-12,
                           weight='bold', fontstyle='italic', fontsize=30)
    # Adjusting the format of the colorbar range
    if exp_form == True:
        xmin_cbar, xmax_cbar = fig.axes[1].get_xlim()
        scale_pow = int(np.log10(xmax_cbar))

        def my_formatter_fun(x, p):
            return "%.2f" % (x / (10 ** scale_pow))
        fig.axes[1].get_xaxis().set_major_formatter(
            ticker.FuncFormatter(my_formatter_fun))
        fig.axes[1].set_xlabel(label, fontsize=26)
        plt.annotate('x $10^{{{0:d}}}$'.format(scale_pow), xy=(0.57, -0.25),
                     fontsize=22, xycoords='axes fraction', color='k', alpha=0.6)
    fig.savefig('all_medial.svg', dpi=150, transparent=True)
    plt.close()

    # Do the dorsal view
    p = Plot(lh, rh, zoom=3, views=[
             'dorsal'], flip=False, brightness=brightness, size=(400, 500))
    p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=False, cmap=cmap,
                color_range=(thresh, upper_thresh))
    if edges == True:
        p.add_layer({'left': lh_regions, 'right': rh_regions}, cmap=mygraymap,
                    as_outline=True, cbar=False)
    fig = p.build(scale=(2, 2))
    fig.savefig('dorsal.svg', dpi=150, transparent=True)
    plt.close()

    # Do the ventral view
    p = Plot(lh, rh, zoom=3, views=[
             'ventral'], flip=True, brightness=brightness, size=(400, 500))

    p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=False, cmap=cmap,
                color_range=(thresh, upper_thresh))
    if edges == True:
        p.add_layer({'left': lh_regions, 'right': rh_regions}, cmap=mygraymap,
                    as_outline=True, cbar=False)

    fig = p.build(scale=(2, 2))
    fig.savefig('ventral.svg', dpi=150, transparent=True)
    plt.close()

    # DO the anterior view
    p = Plot(lh, rh, zoom=3.2, views=[
             'anterior'], flip=True, brightness=brightness, size=(400, 400), mirror_views=True)
    p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=False, cmap=cmap,
                color_range=(thresh, upper_thresh))
    if edges == True:
        p.add_layer({'left': lh_regions, 'right': rh_regions}, cmap=mygraymap,
                    as_outline=True, cbar=False)

    fig = p.build(scale=(2, 2))
    fig.savefig('anterior.svg', dpi=150, transparent=True)
    plt.close()

    # Do the posterior view
    p = Plot(lh, rh, zoom=3.2, views=[
             'posterior'], flip=False, brightness=brightness, size=(400, 400), mirror_views=True)

    p.add_layer({'left': lh_regions, 'right': rh_regions}, cbar=False, cmap=cmap,  # cbar_label='Node Strength',
                color_range=(thresh, upper_thresh))
    if edges == True:
        p.add_layer({'left': lh_regions, 'right': rh_regions}, cmap=mygraymap,
                    as_outline=True, cbar=False)

    fig = p.build(scale=(2, 2))
    fig.savefig('posterior.svg', dpi=150, transparent=True)
    plt.close()


def compose_full_view(remove_files=True, output_name='fig_final', save_svg=False, save_png=False):
    '''
    Assemble the full brain map by considering the 6 views considered
    using the full_view function.
    The function saves the figure as output_name.svg inthe current 
    working directory
    '''
    fig = sg.SVGFigure()
    fig.set_size(('3.75cm', '3cm'))
    fig_bottom_left = sg.fromfile('anterior.svg')
    fig_bottom_right = sg.fromfile('posterior.svg')
    fig_central_bottom = sg.fromfile('ventral.svg')
    fig_central_upper = sg.fromfile('dorsal.svg')
    fig_upper = sg.fromfile('all_lateral.svg')
    fig_bottom = sg.fromfile('all_medial.svg')

    pl_cb = fig_central_bottom.getroot()
    pl_cu = fig_central_upper.getroot()
    pl_bl = fig_bottom_left.getroot()
    pl_br = fig_bottom_right.getroot()
    pl_u = fig_upper.getroot()
    pl_b = fig_bottom.getroot()
    scaling = 0.15
    scaling_central = 0.145
    scaling_medial = 0.172
    scaling_down = 0.135
    pl_cb.rotate(180, x=0, y=0)

    pl_bl.moveto(-6, 68, scale_x=scaling_down, scale_y=scaling_down)
    pl_br.moveto(98, 68, scale_x=scaling_down, scale_y=scaling_down)

    pl_cb.moveto(97.5, 105, scale_x=scaling_central, scale_y=scaling_central)
    pl_cu.moveto(43.5, -8, scale_x=scaling_central, scale_y=scaling_central)

    pl_u.moveto(-52, -16.5, scale_x=scaling_central, scale_y=scaling_central)
    pl_b.moveto(-81.5, 24, scale_x=scaling_medial, scale_y=scaling_medial)

    txt1 = sg.TextElement(0, 5, "L", size=6, color="k")
    txt2 = sg.TextElement(137.5, 5, "R", size=6, color="k")

    fig.append([pl_bl, pl_br, pl_cb, pl_cu, pl_u, pl_b])
    fig.append([txt1, txt2])

    fig.save('{0}.svg'.format(output_name))

    # os.system(
    #     'inkscape --export-png={0}.png {0}.svg -d 600 -b white'.format(output_name))
    os.system(
        'inkscape {0}.svg -o {0}.png  -d 900 -b white'.format(output_name))
    if remove_files == True:
        os.system(
            'rm anterior.svg posterior.svg ventral.svg dorsal.svg all_lateral.svg all_medial.svg')
    if save_svg == False:
        os.system('rm {0}.svg'.format(output_name))

    fig1 = plt.figure(dpi=150)
    im = plt.imread('{0}.png'.format(output_name))
    imgplot = plt.imshow(im)

    if save_png == False:
        os.system('rm {0}.png'.format(output_name))

    ax = plt.gca()
    ax.axis('off')
    return(fig1)
