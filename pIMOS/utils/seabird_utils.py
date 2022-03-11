"""
    Tools for working with Seabird CTD files.

    Some code from UHDAS

"""

import re
import datetime
import numpy as np

from datetime import datetime, timedelta
import xarray as xray
import os

import pdb

##########
# UHDAS code
##########
#import pycurrents.system.logutils as logutils
#L = logutils.getLogger(__file__)

def shift_CTD_yearbase(yearbase, dday):
    """
    Convert TimeQ days (divide by 86400 first to get to days from
    the original seconds) to decimal days with the given *yearbase*.

    There is no guarantee this will work right on other files, but
    it does work for data from the KM in 2012.
    """
    dday_offset = datetime.date(yearbase, 1, 1) - datetime.date(2000, 1, 1)
    return dday - dday_offset.days

def _dedup(names):
    """
    Given a list of names, return a list in which duplicates have
    been munged to make them unique.
    This is needed in case someone has configured SeaSoft to
    produce cnv files with duplicate columns.
    """
    nameset = set(names)
    if len(names) == len(nameset):
        return names

    for i in range(1, len(names)):
        if names[i] in names[:i]:
            #L.warn("Duplicate name in CTD file: %s", names[i])
            names[i] = '_' + names[i]

    return names

def _to_ascii_only(names):
    _names = []
    theta = bytes([0xe9])
    for name in names:
        name_b = name.encode('CP437')
        if theta in name_b:
            head, tail = name_b.split(theta)
            name = (head + b'theta' + tail).decode('CP437')
        name = name.replace('-', '_')
        name = name.replace('/', 'per')
        _names.append(name)
    return _names

class CnvFile(object):
    """
    Seabird CTD cnv file reader.

    Initialize it with the file name;
    access the data as a rectangular array with the 'array' attribute;
    or as a record array with the 'records' attribute.
    """

    IDkeys = ['Cruise ID', 'Station', 'Date', 'Depth(m)']
    NMEAkeys = ['NMEA Longitude', 'NMEA Latitude', 'NMEA UTC (Time)']

    def __init__(self, filename, yearbase=None, edit=True):
        """
        Given a cnv filename, open and parse the file.

        If *edit* is True, the returned data will be truncated
        where the scan count goes backwards.

        This code is derived from Roberto de Almeida's blog post
        of February 7, 2008, here:
        http://swik.net/PyTextile/PyTextile+and+PyDap+Blog/Reading+Seabird+files/b2nm8

        Switched to double precision so that the time fields would work.
        """
        self.filename = filename
        self.yearbase = yearbase # will be filled in below if it is None

        # cnv files normally come with DOS line endings, but sometimes
        # one may need to edit the file, in which case the line endings
        # might be converted to unix endings.  Therefore we handle
        # both cases.

        header = b''
        search_start = 0
        _end = b'\n*END*'
        with open(filename, 'rb') as infile:
            while True:
                chunk = infile.read(1024)
                header += chunk
                i0 = header[search_start:].find(_end)
                if i0 > -1:
                    i0 += search_start
                    if header[i0-1:i0] == b'\r':
                        line_ending = '\r\n'
                        idata = i0 + 8
                    else:
                        line_ending = '\n'
                        idata = i0 + 7
                    header = header[:i0]
                    break
                if len(chunk) < 1024:
                    raise RuntimeError('end of header not found')
                if search_start == 0:
                    search_start = 1024 - len(_end)
                else:
                    search_start += 1024

        self.i0 = i0
        self.idata = idata

        header = header.decode('CP437')  # Ancient DOS character set.

        # Process headers.
        p = re.compile("""
            (?:\*{1,2}|\#)\s    # starts with "* " or "** " or "# "
            (.*?)               # everything until next token
            \s*(?::|=|$)\s*     # ":  " or " = " or EOL
            (.*)                # all the rest, if any
        """, re.VERBOSE)
        lines = header.split(line_ending)
        #attributes = [p.match(line).groups() for line in lines if line.strip()]
        attributes = []
        for line in lines:
            #print line
            if line.strip():
                match  = p.match(line)
                if match is not None:
                    attributes.append(match.groups())

        # Filter out the XML junk.
        attributes = [(k.strip(), v.strip()) for (k, v) in attributes
                                    if not k.strip().startswith('<')]
        self.attributes = dict(attributes)
        self.ordered_keys = [k for (k, v) in attributes]
        self.file_type = self.attributes.get('file_type', 'binary').lower()
        if not self.file_type in ('binary', 'ascii'):
            raise ValueError("file_type is %s; must be ascii or binary" %
                                self.file_type)
        fields = []
        longnames = []
        cols = []
        for a in self.ordered_keys:
            if a.startswith('name '):
                col = int(a.split()[1])
                cols.append(col)        # for a sanity check; may be removed
                field, longname = self.attributes[a].split(':')
                fields.append(field)
                longnames.append(longname)
        self.nfields = int(self.attributes.pop('nquan'))
        # nvalues does not seem to correspond to what is in the file
        #self.nrecords = int(self.attributes.pop('nvalues'))
        self._nrecords = None
        assert len(cols) == self.nfields, "Wrong number of names"
        fields = _dedup(fields)
        longnames =  _dedup(longnames)
        fields = _to_ascii_only(fields)  # Py2 workaround for np.ma problem?
        self._dtype = np.dtype({'names': fields,
                                'formats': ['f8']*self.nfields,
                                'titles': longnames})
        self.names = fields
        self.longnames = longnames
        self._array = None
        self._records = None
        try:
            self.lon = self._decimal_lon(self.attributes['NMEA Longitude'])
            self.lat = self._decimal_lat(self.attributes['NMEA Latitude'])
        except:
            print('No lat/lon stored.')

        # Try to get the date and time from NMEA; sometimes this has
        # date and time, sometimes just the time.
        #nmea_time = self.attributes['NMEA UTC (Time)']
        nmea_time = self.attributes['start_time']
        try:
            tfmt = "%b %d %Y %H:%M:%S"
            self.datetime = datetime.strptime(nmea_time[0:20], tfmt)
        except ValueError:
            # It must have only the time; get the date from Date.
            _date = self.attributes['Date']
            tfmt = "%d-%b-%Y %H:%M:%S"
            self.datetime = datetime.strptime(
                            '%s %s' % (_date, nmea_time), tfmt)

        dt = self.datetime.timetuple()
        self.ymdhms = dt[:6]
        self.yeardaynum = dt[7]
        if self.yearbase is None:
            self.yearbase = self.ymdhms[0]
            y_adjust = 0
        else:
            y_adjust = (datetime.date(self.ymdhms[0], 1, 1)
                        -datetime.date(self.yearbase, 1, 1)).days
        h,m,s = self.ymdhms[3:]
        self.dday_start = (self.yeardaynum - 1) + h/24.0 + m/1440.0 + s/86400.0
        self.dday_start += y_adjust
        self._NMEA_dday = None
        self._scan_dday = None
        self._dday = None
        self._pressure = None

        self.edit = edit


    @property
    def scan_dday(self):
        """
        dday based on scan time ('timeS') instead of NMEA ('timeQ')
        """
        if self._dday is None:
            d = self.records['timeS']/86400
            #if d.count() != d.size:
            #    L.warn("Some scan times were masked and are now NaN;"
            #            " scan_dday is unreliable.")
            self._scan_dday = self.dday_start + d.filled(np.nan)
        return self._scan_dday


    @property
    def NMEA_dday(self):
        """
        dday based on NMEA times.
        """

        if self._NMEA_dday is None:
            if 'timeQ' not in self.names:
                d = np.empty((len(self.records),), dtype=np.float)
                d.fill(np.nan)
                #L.warn("This file has no 'timeQ'.")
                self._NMEA_dday = d
            else:
                d = self.records['timeQ'] / 86400.0
                #if d.count() != d.size:
                #    L.warn("Some NMEA times were masked and are now NaN;"
                #            " NMEA_dday is unreliable.")
                d = d.filled(np.nan)
                self._NMEA_dday = shift_CTD_yearbase(self.yearbase, d)

        return self._NMEA_dday

    @property
    def dday(self):
        """
        Attempt to find or generate a monotonic dday vector.
        """
        if self._dday is None:
            d = self.NMEA_dday
            if np.isnan(d).any() or (np.diff(d) <= 0).any():
                #L.info("Rejected NMEA_dday")
                d = self.scan_dday
                kind, interval = self.attributes['interval'].split()
                if (kind != "decibars:" and
                        (np.isnan(d).any() or (np.diff(d) <= 0).any())):
                    #L.info("Rejected scan_dday, using scans")
                    d = self.records['scan'] * float(interval)
                    d = self.dday_start + d / 86400.0
            self._dday = d
        return self._dday

    # The following assumes pressure is always prDM; not sure if it is true.
    @property
    def pressure(self):
        if self._pressure is None:
            self._pressure = self.records['prDM']
        return self._pressure

    def __str__(self):
        sl = ["filename: %s" % self.filename]
        for key in self.IDkeys:
            try:
                sl.append("%s:  %s" % (key, self.attributes[key]))
            except KeyError:
                pass
        for key in self.NMEAkeys:
            sl.append("%s:  %s" % (key, self.attributes[key]))
        sl.append("%d records of %d fields each" %
                        (self.nrecords, self.nfields))
        sl.append("field names:")
        sl.append("\n".join(self.names))
        return "\n".join(sl)

    @staticmethod
    def _decimal_lon(xstr):
        deg, min, EW = xstr.split()
        d = float(deg) + float(min) / 60.0
        if EW == 'W':
            d = -d
        return d

    @staticmethod
    def _decimal_lat(ystr):
        deg, min, NS = ystr.split()
        d = float(deg) + float(min) / 60.0
        if NS == 'S':
            d = -d
        return d

    def _get_array(self):
        if self._array is None:
            with open(self.filename, 'rb') as infile:
                infile.seek(self.idata)
                _databytes = infile.read()

            if self.file_type == 'binary':
                data = np.fromstring(_databytes, '<f4')
                data.shape = (-1, self.nfields)
                if 'timeQ' in self.names:
                    i_timeQ = self.names.index('timeQ')
                    t = data[:, i_timeQ].view(np.uint32)
                    data = data.astype(float)
                    data[:, i_timeQ] = t
            else:
                strdata = _databytes.decode('ascii')
                data = np.fromstring(strdata,
                                     sep=' ', dtype='f8')
                data.shape = (-1, self.nfields)

            if self.edit:
                scan = data[:,0]
                dscan = np.diff(scan)
                glitch = np.nonzero(dscan <=0)[0]
                if glitch.size > 0:
                    #L.warn("Found %s nonmonotonic scan points.", glitch.size)
                    i1 = glitch[0]
                    n0 = data.shape[0]
                    data = data[:i1]
                    #L.warn("Truncated %s at index %s of %s.",
                    #                    self.filename, i1, n0)
            data = np.ma.array(data)
            bad_value = float(self.attributes['bad_flag'])
            cond = np.ma.abs(data - bad_value) < 0.01 * abs(bad_value)
            data[cond] = np.ma.masked
            # This should work if bad_value is consistently something
            # like the present tiny negative number.
            self._array = data
        
        return self._array

    array = property(_get_array)

    def _get_records(self):
        if self._records is None:
            a = self._get_array()
            self._records = a.view(dtype=self._dtype)
            self._records.shape = (a.shape[0],)
        return self._records

    records = property(_get_records)

    def _get_nrecords(self):
        if self._nrecords is None:
            self._nrecords = self.array.shape[0]
        return self._nrecords

    nrecords = property(_get_nrecords)
