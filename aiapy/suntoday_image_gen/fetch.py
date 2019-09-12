import requests
import datetime
import dateutil
import pandas as pd
from io import StringIO


def getKeys(dataseries, useJSOC2=False):
    """Get valid key names."""
    params = {
        'ds': dataseries + '[$]',
        'a': 1
    }
    text = _getText(_getURL(useJSOC2), params)
    return text.split('\n')[0].split('\t')


def getSegments(dataseries, useJSOC2=False):
    """Get valid segment names."""
    params = {
        'ds': dataseries + '[$]',
        'A': 1
    }
    text = _getText(_getURL(useJSOC2), params)
    return text.split('\n')[0].split('\t')


def fetch(dataseries, start=None, end_or_span=None, cadence=None,
          wavelengths=None, segments=None, keys=None, query=False,
          useJSOC2=False, headers=False, df=False, parse_dates=False):
    """Fetch file locations.

    Parameters
    ----------
    dataseries : 'aia' or 'hmi.sharp' or 'hmi.B' or str
        The dataseries to select from JSOC.
    start : datetime or str, optional
        Starting datetime or a string to parse as a datetime.
    end_or_span : datetime or timedelta or int or str, optional
        Either the ending date or the time span over which to search the data.
            If datetime, treat as ending time.
            If timedelta, use as time span.
            If int, treat as number of seconds.
            If str, directly place in URL as span parameter.
    cadence : timedelta or int or str, optional
        The cadence at which to search the data. Same treatment as end_or_span,
        not including datetime.
    wavelengths : int or sequence of ints, optional
        Wavelength(s) to search.
    segments : str or sequence of strs or 'all', optional
        Return the file paths for this/these segment(s).
    keys : str or sequence of strs or 'all', optional
        Return the values for this/these key(s).
    query : boolean, optional
        Return the full query for each row.
    useJSOC2 : boolean, optional
        Use the JSOC2 url instead of JSOC.
    headers : boolean, optional
        Whether to return column names preceding the rows in the returned list.
        Has no effect if df=True.
    df : boolean, optional
        Whether to return a pandas.DataFrame or list of strings.
    parse_dates : list of strs, optional
        List of columns to parse as datetimes before returning dataframe. Has
        no effect if df=False.

    Returns
    -------
    A pandas.DataFrame object or a list of returned data rows.
    """
    if start:
        if isinstance(start, str):
            start = dateutil.parser.parse(start)
        if dataseries[:3] == 'aia':
            dataseries += '[{:%Y-%m-%dT%H:%M:%S}Z'.format(start)
        elif dataseries[:3] == 'hmi':
            dataseries += '[{:%Y.%m.%d_%H:%M:%S}_TAI'.format(start)
        else:
            dataseries += '[{}'.format(start)

        if end_or_span:
            dataseries += '/' + _convertToSeconds(end_or_span, start)
            if cadence:
                dataseries += '@' + _convertToSeconds(cadence)
        dataseries += ']'

    if wavelengths:
        if not isinstance(wavelengths, int):
            wavelengths = ','.join(str(c) for c in wavelengths)
        dataseries += '[{}]'.format(wavelengths)

    params = {'ds': dataseries}

    if segments:
        params['P'] = 1  # Show full file path
        if segments == 'all':
            params['A'] = 1
        else:
            if not isinstance(segments, str):
                segments = ','.join(segments)
            params['seg'] = segments
    if keys:
        if keys == 'all':
            params['a'] = 1
        else:
            if not isinstance(keys, str):
                keys = ','.join(keys)
            params['key'] = keys
    if query:
        params['i'] = 1

    if df:
        text = _getText(_getURL(useJSOC2), params)
        if text:
            return pd.read_csv(StringIO(text),
                               sep='\t', parse_dates=parse_dates)
        else:
            return pd.DataFrame()
    else:
        if not headers:
            params['q'] = 1
        text = _getText(_getURL(useJSOC2), params)

        if '\t' in text:
            return [line.split('\t') for line in text.split('\n')[:-1]]
        else:
            return text.split('\n')[:-1]


def _getText(URL, params):
    """Get the text of JSOC Lookdata, but raise a RuntimeError if an error
    message is present.
    """
    r = requests.get(URL, params=params, headers={'connection': 'close'})
    if 'Error' in r.text:
        raise RuntimeError(
            'JSOC Lookdata failed.\nURL: {}\nOutput:\n\n{}'.format(
                r.url, r.text))
    else:
        return r.text


def _convertToSeconds(timeSpan, start=None):
    """Convert timeSpan to number of seconds."""
    if isinstance(timeSpan, str):
        return timeSpan
    else:
        if isinstance(timeSpan, datetime.datetime):
            seconds = (timeSpan - start).total_seconds()
        elif isinstance(timeSpan, datetime.timedelta):
            seconds = timeSpan.total_seconds()
        else:
            seconds = timeSpan
        return str(seconds) + 's'


def _getURL(useJSOC2):
    """Return either the original JSOC url or the JSOC2 URL."""
    two = '2' if useJSOC2 else ''
    return 'http://jsoc{}.stanford.edu/cgi-bin/ajax/show_info'.format(two)
