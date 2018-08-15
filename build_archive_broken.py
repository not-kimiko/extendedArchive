from __future__ import absolute_import, division, print_function
import os
import re
import yaml
import argparse
import xml.etree.cElementTree as ElementTree
import numpy as np
from astropy.io.fits import FITS_rec, Column, BinTableHDU
from astropy.coordinates import SkyCoord


def path_to_xmlpath(path):
    if path is None:
        return path

    return re.sub(r'\$([a-zA-Z\_]+)', r'$(\1)', path)


def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir


def prettify_xml(elem):
    """Return a pretty-printed XML string for the Element.
    """
    from xml.dom import minidom
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


def isstr(s):
    """String instance testing method that works under both Python 2.X
    and 3.X.  Returns true if the input is a string."""

    try:
        return isinstance(s, basestring)
    except NameError:
        return isinstance(s, str)


def to_xml(xmlfile, name, src_dict):

    params = src_dict['Spectral_Parameters']
    spatial_type = src_dict['Spatial_Function']
    spectral_type = src_dict['Spectral_Function']
    defaults = {
        'Prefactor': {'scale': 1, 'name': 'Prefactor', 'free': 0, 'min': 1, 'max': 1},
        'RA': {'scale': 1, 'name': 'RA', 'free': 0, 'min': 0, 'max': 360},
        'DEC': {'scale': 1, 'name': 'DEC', 'free': 0, 'min': -90, 'max': 90},
        'Radius': {'scale': 1, 'name': 'Radius', 'free': 0, 'min': 0, 'max': 10},
        'Sigma': {'scale': 1, 'name': 'Sigma', 'free': 0, 'min': 0, 'max': 10},
    }

    if spatial_type == 'RadialDisk':
        params_spat = {'RA': defaults['RA'].copy(), 'DEC': defaults['DEC'].copy(),
                       'Radius': defaults['Radius'].copy()}
        params_spat['Radius']['value'] = np.sqrt(src_dict['Model_SemiMajor'] *
                                                 src_dict['Model_SemiMinor'])
        params_spat['RA']['value'] = src_dict['RAJ2000']
        params_spat['DEC']['value'] = src_dict['DEJ2000']
    elif spatial_type == 'RadialGaussian':
        params_spat = {'RA': defaults['RA'].copy(), 'DEC': defaults['DEC'].copy(),
                       'Sigma': defaults['Sigma'].copy()}
        params_spat['Sigma']['value'] = np.sqrt(src_dict['Model_SemiMajor'] *
                                                src_dict['Model_SemiMinor'])
        params_spat['Sigma']['value'] /= 1.5095921854516636
        params_spat['RA']['value'] = src_dict['RAJ2000']
        params_spat['DEC']['value'] = src_dict['DEJ2000']
    else:
        params_spat = { 'Prefactor' : defaults['Prefactor'] }

    root = ElementTree.Element('source_library')
    root.set('title', 'source_library')

    source_element = create_xml_element(root, 'source',
                                        dict(name = name,
                                             Photon_Flux = src_dict['Photon_Flux'],
                                             Energy_Flux = src_dict['Energy_Flux'],
                                             type='DiffuseSource'))

    #filename = utils.path_to_xmlpath(self.filefunction)
    spec = create_xml_element(source_element, 'spectrum',
                              dict(type=src_dict['Spectral_Function']))

    attrib = dict(type=spatial_type)
    if spatial_type == 'SpatialMap':
        attrib['file'] = path_to_xmlpath(src_dict['Spatial_Filename'])
        attrib['map_based_integral'] = 'true'
    
    spat = create_xml_element(source_element, 'spatialModel', attrib)

    for k, v in params.items():
        create_xml_element(spec, 'parameter', v)

    for k, v in params_spat.items():
        create_xml_element(spat, 'parameter', v)

    output_file = open(xmlfile, 'w')
    output_file.write(prettify_xml(root))


def create_xml_element(root, name, attrib):
    el = ElementTree.SubElement(root, name)
    for k, v in attrib.items():

        if isinstance(v, bool):
            el.set(k, str(int(v)))
        elif isstr(v):
            el.set(k, v)
        elif np.isfinite(v):
            el.set(k, str(v))

    return el


def set_coordinates(src_dict):

    has_cel = 'RAJ2000' in src_dict and 'DEJ2000' in src_dict
    has_gal = 'GLON' in src_dict and 'GLAT' in src_dict

    if not has_cel and not has_gal:
        raise Exception()

    if not has_gal:
        skydir = SkyCoord(src_dict['RAJ2000'], src_dict['RAJ2000'],
                          unit='deg', frame='icrs')
        skydir = skydir.transform_to('galactic')
        src_dict['GLAT'] = skydir.b.deg
        src_dict['GLON'] = skydir.l.deg


def build_column_array(column_name,sources,npar_max):
    spec_par_list = [
        'Spectral_Param_Name',
        'Spectral_Param_Value',
        'Spectral_Param_Error',
        'Spectral_Param_Scale'
        ]

    if column_name in spec_par_list:
    
    	for k, v in sources.items():
    	
            params = v.get('Spectral_Parameters')
            npars = len(params.keys())
    
        if column_name == 'Spectral_Param_Name':
            collist = [str(params[i]['name']) for i in params.keys()] + [''] * (npar_max - npars)
            
        elif column_name == 'Spectral_Param_Value':
	    collist = [float(params[i]['value']) for i in params.keys()] + [np.nan] * (npar_max - npars)
			
	elif column_name == 'Spectral_Param_Scale':
	    collist = [float(params[i]['scale']) for i in params.keys()] + [np.nan] * (npar_max - npars)
			
	else:
            collist = [float(params[i]['error']) if 'error' in params[i] else np.nan for i in params.keys()] \
				+ [np.nan] * (npar_max - npars)

    else: 
	collist = [sources[k].get(column_name, None) for k in sources.keys()]

    return np.asarray(collist)


def main():

    usage = "usage: %(prog)s [archive file]"
    description = "Build the extended archive from the master archive YAML file."
    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('--outname', default=None, required=True)

    parser.add_argument('--vernum', default=0, required=True)

    parser.add_argument('masterfile',
                        help='Extended archive master YAML file.')

    args = parser.parse_args()
    
    npar_max = 5
    sources = yaml.load(open(args.masterfile))
    cols = [Column(name='Source_Name', format='18A'),
            Column(name='RAJ2000', format='E', unit='deg', disp='F8.4'),
            Column(name='DEJ2000', format='E', unit='deg', disp='F8.4'),
            Column(name='GLON', format='E', unit='deg', disp='F8.4'),
            Column(name='GLAT', format='E', unit='deg', disp='F8.4'),
            Column(name='Photon_Flux', format='E', unit='ph cm-2 s-1', disp='E8.2'),
            Column(name='Energy_Flux', format='E', unit='erg cm-2 s-1', disp='E8.2'),
            Column(name='Model_Form', format='12A'),
            Column(name='Model_SemiMajor', format='E', unit='deg', disp='E7.3'),
            Column(name='Model_SemiMinor', format='E', unit='deg', disp='E7.3'),
            Column(name='Model_PosAng',  format='E', unit='deg', disp='E6.1'),
            Column(name='Spatial_Function', format='15A'),
            Column(name='Spatial_Filename', format='50A'),
            Column(name='Spectral_Function', format='12A'),
            Column(name='Spectral_Filename', format='40A'),
            Column(name='Name_1FGL', format='18A'),
            Column(name='Name_2FGL', format='18A'),
            Column(name='Name_3FGL', format='18A'),
#            Column(name='Spectral_Param_Name', format='45A9'),
#            Column(name='Spectral_Param_Value',
#                    format='E', dim=str(npar_max), disp='E9.4'),
#            Column(name='Spectral_Param_Error',
#                    format='E', dim=str(npar_max), disp='E9.4'),
#            Column(name='Spectral_Param_Scale',
#                    format='E', dim=str(npar_max)),
            ]


    for c in cols: 
        c.array = build_column_array(c.name, sources, npar_max)
 
    
    record = FITS_rec.from_columns(cols)
       
       
    outdir = args.outname + "_v" + args.vernum
    mkdir(outdir)

    fitsname = "LAT_extended_sources_v" + args.vernum + ".fits"
    output = BinTableHDU(record)
    output.writeto(os.path.join(outdir,
                           fitsname), overwrite=True)


    xmldir = os.path.join(outdir, 'XML')
    mkdir(xmldir)
    
    for k, v in sources.items():
        xmlpath = os.path.join(xmldir, v['Source_Name'].replace(' ', '') + '.xml')
        to_xml(xmlpath, v['Source_Name'], v)



    # TODO: Generate region file


if __name__ == "__main__":
    main()
