#!/usr/bin/env python
# encoding: utf-8
"""
core.py

Created by José Sánchez-Gallego on 30 Sep 2013.
Copyright (c) 2013. All rights reserved.

Classes and functions to analyse the PDR results.
"""

__all__ = ['ConfigureAnaliser']


import ConfigParser
from warnings import warn, simplefilter, resetwarnings
from ..exceptions import AnalysePDRDataError


class ConfigureAnaliser(object):

    def __init__(self, configFile):

        if not isinstance(configFile, str):
            raise TypeError('configFile is not of type str')

        self.configFile = configFile
        self.imageData = {}

        self._parseConfigFile()

        if len(self.imageData.keys()) == 0:
            raise AnalysePDRDataError('No images to analyse')

        return

    def _parseConfigFile(self):

        configFile = self.configFile

        configData = ConfigParser.SafeConfigParser()
        configData.read(configFile)
        configData.optionxform = str

        fileFields = configData.options('Images')

        for field in fileFields:

            if 'image' in field.lower():
                image = configData.get('Images', field)
                imageId = field[len('image'):]

                conversionField = 'conversion' + imageId
                if configData.has_option('Images', conversionField):
                    try:
                        conversionFactor = configData.getfloat('Images', conversionField)
                    except:
                        raise TypeError('conversion%s is not a float' % imageId)
                else:
                    raise ConfigParser.NoOptionError('Field %s cannot be found' % conversionField)

                unitField = 'unit%s' % imageId
                if configData.has_option('Images', unitField):
                    unit = configData.get('Images', unitField)
                else:
                    unit = None

                isImageOK = self.checkImage(image)
                if isImageOK is True:
                    self.imageData[image] = {'conversionFactor': conversionFactor,
                                             'unit': unit.strip()}

        return

    def checkImage(self, image):
        """
        Checks if an image is a fits file with a header
        and valid WCS information.
        """

        from astropy.io import fits
        from astropy import wcs

        simplefilter('always')
        try:
            hdu = fits.open(image)
            header = hdu[0].header
            imageWCS = wcs.WCS(header)
        except:
            warn('Image %s cannot be read or has no astrometry. It wont be used' % image)
            return False
        resetwarnings()
        del hdu, header, imageWCS

        return True
