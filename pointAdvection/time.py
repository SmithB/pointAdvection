#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
time.py
Written by Tyler Sutterley (05/2022)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/

UPDATE HISTORY:
    Updated 05/2022: changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: date parser for cases when only a calendar date
    Written 07/2020
"""
import os
import re
import copy
import logging
import warnings
import datetime
import numpy as np
import dateutil.parser

#-- PURPOSE: parse a date string into epoch and units scale
def parse_date_string(date_string):
    """
    parse a date string of the form

    - time-units since ``yyyy-mm-dd hh:mm:ss``
    - ``yyyy-mm-dd hh:mm:ss`` for exact calendar dates

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss

    Returns
    -------
    epoch: list
        epoch of delta time
    conversion_factor: float
        multiplication factor to convert to seconds
    """
    #-- try parsing the original date string as a date
    try:
        epoch = dateutil.parser.parse(date_string)
    except ValueError:
        pass
    else:
        #-- return the epoch (as list)
        return (datetime_to_list(epoch),0.0)
    #-- split the date string into units and epoch
    units,epoch = split_date_string(date_string)
    conversion_factors = {'microseconds': 1e-6,'microsecond': 1e-6,
        'microsec': 1e-6,'microsecs': 1e-6,
        'milliseconds': 1e-3,'millisecond': 1e-3,'millisec': 1e-3,
        'millisecs': 1e-3,'msec': 1e-3,'msecs': 1e-3,'ms': 1e-3,
        'seconds': 1.0,'second': 1.0,'sec': 1.0,'secs': 1.0,'s': 1.0,
        'minutes': 60.0,'minute': 60.0,'min': 60.0,'mins': 60.0,
        'hours': 3600.0,'hour': 3600.0,'hr': 3600.0,
        'hrs': 3600.0,'h': 3600.0,
        'day': 86400.0,'days': 86400.0,'d': 86400.0}
    if units not in conversion_factors.keys():
        raise ValueError('Invalid units: {0}'.format(units))
    #-- return the epoch (as list) and the time unit conversion factors
    return (datetime_to_list(epoch),conversion_factors[units])

#-- PURPOSE: split a date string into units and epoch
def split_date_string(date_string):
    """
    split a date string into units and epoch

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss
    """
    try:
        units,_,epoch = date_string.split(None,2)
    except ValueError:
        raise ValueError('Invalid format: {0}'.format(date_string))
    else:
        return (units.lower(),dateutil.parser.parse(epoch))

#-- PURPOSE: convert a datetime object into a list
def datetime_to_list(date):
    """
    convert a datetime object into a list

    Parameters
    ----------
    date: datetime object

    Returns
    -------
    date: list
        [year,month,day,hour,minute,second]
    """
    return [date.year,date.month,date.day,date.hour,date.minute,date.second]

#-- PURPOSE: gets the number of days per month for a given year
def calendar_days(year):
    """
    Calculates the number of days per month for a given year

    Parameters
    ----------
    year: int or float
        calendar year

    Returns
    -------
    dpm: list
        number of days for each month
    """
    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31],dtype=np.float64)
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float64)
    #-- Rules in the Gregorian calendar for a year to be a leap year:
    #-- divisible by 4, but not by 100 unless divisible by 400
    #-- True length of the year is about 365.2422 days
    #-- Adding a leap day every four years ==> average 365.25
    #-- Subtracting a leap year every 100 years ==> average 365.24
    #-- Adding a leap year back every 400 years ==> average 365.2425
    #-- Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (year % 4)
    m100 = (year % 100)
    m400 = (year % 400)
    m4000 = (year % 4000)
    #-- find indices for standard years and leap years using criteria
    if ((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0)):
        return dpm_leap
    elif ((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0)):
        return dpm_stnd

#-- PURPOSE: convert a numpy datetime array to delta times from the UNIX epoch
def convert_datetime(date, epoch=(1970,1,1,0,0,0)):
    """
    Convert a numpy datetime array to seconds since ``epoch``

    Parameters
    ----------
    date: obj
        numpy datetime array
    epoch: tuple, default (1970,1,1,0,0,0)
        epoch for output delta_time

    Returns
    -------
    delta_time: float
        seconds since epoch
    """
    epoch = datetime.datetime(*epoch)
    return (date - np.datetime64(epoch)) / np.timedelta64(1, 's')

#-- PURPOSE: convert times from seconds since epoch1 to time since epoch2
def convert_delta_time(delta_time, epoch1=None, epoch2=None, scale=1.0):
    """
    Convert delta time from seconds since ``epoch1`` to time since ``epoch2``

    Parameters
    ----------
    delta_time: float
        seconds since epoch1
    epoch1: tuple or NoneType, default None
        epoch for input delta_time
    epoch2: tuple or NoneType, default None
        epoch for output delta_time
    scale: float, default 1.0
        scaling factor for converting time to output units
    """
    epoch1 = datetime.datetime(*epoch1)
    epoch2 = datetime.datetime(*epoch2)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    #-- subtract difference in time and rescale to output units
    return scale*(delta_time - delta_time_epochs)

#-- PURPOSE: calculate the delta time from calendar date
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
def convert_calendar_dates(year, month, day, hour=0.0, minute=0.0, second=0.0,
    epoch=(1992,1,1,0,0,0), scale=1.0):
    """
    Calculate the time in time units since ``epoch`` from calendar dates

    Parameters
    ----------
    year: float
        calendar year
    month: float
        month of the year
    day: float
        day of the month
    hour: float, default 0.0
        hour of the day
    minute: float, default 0.0
        minute of the hour
    second: float, default 0.0
        second of the minute
    epoch: tuple, default (1992,1,1,0,0,0)
        epoch for output delta_time
    scale: float, default 1.0
        scaling factor for converting time to output units

    Returns
    -------
    delta_time: float
        days since epoch
    """
    #-- calculate date in Modified Julian Days (MJD) from calendar date
    #-- MJD: days since November 17, 1858 (1858-11-17T00:00:00)
    MJD = 367.0*year - np.floor(7.0*(year + np.floor((month+9.0)/12.0))/4.0) - \
        np.floor(3.0*(np.floor((year + (month - 9.0)/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0*month/9.0) + day + hour/24.0 + minute/1440.0 + \
        second/86400.0 + 1721028.5 - 2400000.5
    epoch1 = datetime.datetime(1858,11,17,0,0,0)
    epoch2 = datetime.datetime(*epoch)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    #-- return the date in units since epoch
    return scale*np.array(MJD - delta_time_epochs/86400.0,dtype=np.float64)

#-- PURPOSE: Converts from calendar dates into decimal years
def convert_calendar_decimal(year, month, day=None, hour=None, minute=None,
    second=None, DofY=None):
    """
    Converts from calendar date into decimal years taking into
    account leap years

    Parameters
    ----------
    year: float
        calendar year
    month: float
        calendar month
    day: float or NoneType, default None
        day of the month
    hour: float or NoneType, default None
        hour of the day
    minute: float or NoneType, default None
        minute of the hour
    second: float or NoneType, default None
        second of the minute
    DofY: float or NoneType, default None
        day of the year (January 1 = 1)

    Returns
    -------
    t_date: float
        date in decimal-year format

    References
    ----------
    .. [1] Dershowitz, N. and E.M. Reingold. 2008.
        Calendrical Calculations.
        Cambridge: Cambridge University Press.
    """

    #-- number of dates
    n_dates = len(np.atleast_1d(year))

    #-- create arrays for calendar date variables
    cal_date = {}
    cal_date['year'] = np.zeros((n_dates))
    cal_date['month'] = np.zeros((n_dates))
    cal_date['day'] = np.zeros((n_dates))
    cal_date['hour'] = np.zeros((n_dates))
    cal_date['minute'] = np.zeros((n_dates))
    cal_date['second'] = np.zeros((n_dates))
    #-- day of the year
    cal_date['DofY'] = np.zeros((n_dates))

    #-- remove singleton dimensions and use year and month
    cal_date['year'][:] = np.squeeze(year)
    cal_date['month'][:] = np.squeeze(month)

    #-- create output date variable
    t_date = np.zeros((n_dates))

    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap=np.array([31,29,31,30,31,30,31,31,30,31,30,31], dtype=np.float64)
    dpm_stnd=np.array([31,28,31,30,31,30,31,31,30,31,30,31], dtype=np.float64)

    #-- Rules in the Gregorian calendar for a year to be a leap year:
    #-- divisible by 4, but not by 100 unless divisible by 400
    #-- True length of the year is about 365.2422 days
    #-- Adding a leap day every four years ==> average 365.25
    #-- Subtracting a leap year every 100 years ==> average 365.24
    #-- Adding a leap year back every 400 years ==> average 365.2425
    #-- Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (cal_date['year'] % 4)
    m100 = (cal_date['year'] % 100)
    m400 = (cal_date['year'] % 400)
    m4000 = (cal_date['year'] % 4000)
    #-- find indices for standard years and leap years using criteria
    leap, = np.nonzero((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0))
    stnd, = np.nonzero((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0))

    #-- calculate the day of the year
    if DofY is not None:
        #-- if entered directly as an input
        #-- remove 1 so day 1 (Jan 1st) = 0.0 in decimal format
        cal_date['DofY'][:] = np.squeeze(DofY)-1
    else:
        #-- use calendar month and day of the month to calculate day of the year
        #-- month minus 1: January = 0, February = 1, etc (indice of month)
        #-- in decimal form: January = 0.0
        month_m1 = np.array(cal_date['month'],dtype=int) - 1

        #-- day of month
        if day is not None:
            #-- remove 1 so 1st day of month = 0.0 in decimal format
            cal_date['day'][:] = np.squeeze(day)-1.0
        else:
            #-- if not entering days as an input
            #-- will use the mid-month value
            cal_date['day'][leap] = dpm_leap[month_m1[leap]]/2.0
            cal_date['day'][stnd] = dpm_stnd[month_m1[stnd]]/2.0

        #-- create matrix with the lower half = 1
        #-- this matrix will be used in a matrix multiplication
        #-- to calculate the total number of days for prior months
        #-- the -1 will make the diagonal == 0
        #-- i.e. first row == all zeros and the
        #-- last row == ones for all but the last element
        mon_mat=np.tri(12,12,-1)
        #-- using a dot product to calculate total number of days
        #-- for the months before the input date
        #-- basically is sum(i*dpm)
        #-- where i is 1 for all months < the month of interest
        #-- and i is 0 for all months >= the month of interest
        #-- month of interest is zero as the exact days will be
        #-- used to calculate the date

        #-- calculate the day of the year for leap and standard
        #-- use total days of all months before date
        #-- and add number of days before date in month
        cal_date['DofY'][stnd] = cal_date['day'][stnd] + \
            np.dot(mon_mat[month_m1[stnd],:],dpm_stnd)
        cal_date['DofY'][leap] = cal_date['day'][leap] + \
            np.dot(mon_mat[month_m1[leap],:],dpm_leap)

    #-- hour of day (else is zero)
    if hour is not None:
        cal_date['hour'][:] = np.squeeze(hour)

    #-- minute of hour (else is zero)
    if minute is not None:
        cal_date['minute'][:] = np.squeeze(minute)

    #-- second in minute (else is zero)
    if second is not None:
        cal_date['second'][:] = np.squeeze(second)

    #-- calculate decimal date
    #-- convert hours, minutes and seconds into days
    #-- convert calculated fractional days into decimal fractions of the year
    #-- Leap years
    t_date[leap] = cal_date['year'][leap] + \
        (cal_date['DofY'][leap] + cal_date['hour'][leap]/24. + \
        cal_date['minute'][leap]/1440. + \
        cal_date['second'][leap]/86400.)/np.sum(dpm_leap)
    #-- Standard years
    t_date[stnd] = cal_date['year'][stnd] + \
        (cal_date['DofY'][stnd] + cal_date['hour'][stnd]/24. + \
        cal_date['minute'][stnd]/1440. + \
        cal_date['second'][stnd]/86400.)/np.sum(dpm_stnd)

    return t_date

#-- PURPOSE: Converts from Julian day to calendar date and time
def convert_julian(JD, **kwargs):
    """
    Converts from Julian day to calendar date and time

    Parameters
    ----------
    JD: float
        Julian Day (days since 01-01-4713 BCE at 12:00:00)
    astype: str or NoneType, default None
        convert output to variable type
    format: str, default 'dict'
        format of output variables

            - ``'dict'``: dictionary with variable keys
            - ``'tuple'``: tuple in most-to-least-significant order
            - ``'zip'``: aggregated variable sets

    Returns
    -------
    year: float
        calendar year
    month: float
        calendar month
    day: float
        day of the month
    hour: float
        hour of the day
    minute: float
        minute of the hour
    second: float
        second of the minute

    References
    ----------
    .. [1] "Numerical Recipes in C", by William H. Press,
        Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
        Cambridge University Press, 1988 (second printing).
    .. [2] Hatcher, D. A., "Simple Formulae for Julian Day Numbers and
        Calendar Dates", Quarterly Journal of the Royal Astronomical
        Society, 25(1), 1984.
    """
    #-- set default keyword arguments
    kwargs.setdefault('astype', None)
    kwargs.setdefault('format', 'dict')
    #-- raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(ASTYPE='astype',FORMAT='format')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn("""Deprecated keyword argument {0}.
                Changed to '{1}'""".format(old,new),
                DeprecationWarning)
            #-- set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    #-- convert to array if only a single value was imported
    if (np.ndim(JD) == 0):
        JD = np.atleast_1d(JD)
        SINGLE_VALUE = True
    else:
        SINGLE_VALUE = False

    JDO = np.floor(JD + 0.5)
    C = np.zeros_like(JD)
    #-- calculate C for dates before and after the switch to Gregorian
    IGREG = 2299161.0
    ind1, = np.nonzero(JDO < IGREG)
    C[ind1] = JDO[ind1] + 1524.0
    ind2, = np.nonzero(JDO >= IGREG)
    B = np.floor((JDO[ind2] - 1867216.25)/36524.25)
    C[ind2] = JDO[ind2] + B - np.floor(B/4.0) + 1525.0
    #-- calculate coefficients for date conversion
    D = np.floor((C - 122.1)/365.25)
    E = np.floor((365.0 * D) + np.floor(D/4.0))
    F = np.floor((C - E)/30.6001)
    #-- calculate day, month, year and hour
    DAY = np.floor(C - E + 0.5) - np.floor(30.6001*F)
    MONTH = F - 1.0 - 12.0*np.floor(F/14.0)
    YEAR = D - 4715.0 - np.floor((7.0+MONTH)/10.0)
    HOUR = np.floor(24.0*(JD + 0.5 - JDO))
    #-- calculate minute and second
    G = (JD + 0.5 - JDO) - HOUR/24.0
    MINUTE = np.floor(G*1440.0)
    SECOND = (G - MINUTE/1440.0) * 86400.0

    #-- convert all variables to output type (from float)
    if kwargs['astype'] is not None:
        YEAR = YEAR.astype(kwargs['astype'])
        MONTH = MONTH.astype(kwargs['astype'])
        DAY = DAY.astype(kwargs['astype'])
        HOUR = HOUR.astype(kwargs['astype'])
        MINUTE = MINUTE.astype(kwargs['astype'])
        SECOND = SECOND.astype(kwargs['astype'])

    #-- if only a single value was imported initially: remove singleton dims
    if SINGLE_VALUE:
        YEAR = YEAR.item(0)
        MONTH = MONTH.item(0)
        DAY = DAY.item(0)
        HOUR = HOUR.item(0)
        MINUTE = MINUTE.item(0)
        SECOND = SECOND.item(0)

    #-- return date variables in output format (default python dictionary)
    if (kwargs['format'] == 'dict'):
        return dict(year=YEAR, month=MONTH, day=DAY,
            hour=HOUR, minute=MINUTE, second=SECOND)
    elif (kwargs['format'] == 'tuple'):
        return (YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
    elif (kwargs['format'] == 'zip'):
        return zip(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)