# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SiImporter
                                 A QGIS plugin
 Si maps import/export
                             -------------------
        begin                : 2016-06-21
        copyright            : (C) 2016 by Olof Astrand
        email                : olof.astrand@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load SiImporter class from file SiImporter.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .si_importer import SiImporter
    return SiImporter(iface)
