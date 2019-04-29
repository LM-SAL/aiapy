>> Module for the Heliophysics Events Knowledgebase <<
>> Contributors: Nicholas Padmanabhan, Mark Cheung <<

File: HEK.py
Class: HER_Event

The purpose of this class is to generate xml files for the iSolSearch database at lmsal.com/isolsearch. The user creates a HER_Event object, edits its data dictionary with specific values pertaining to the solar event, and uses its export method to save an xml file.

>> Methods <<
1. __init__(id) - takes an event code (see below) and instantiates a HER_Event object with the appropiate "required" and "optional" dictionaries

2. exportEvent(filename = None, suffix = "") - creates an xml file containing all the inputted data

3. __str__() - prints a string representation of an HER_Event object; called with Python's built-in "print" method

>> Example usage (from http://lmsal.com/hek/api.html) <<
Note: If one or more required fields have not been filled out by the user, the xml will not export, and a message will be printed saying which fields are missing.

Note: print(event) can be used to view the current state of the HER_Event dictionary, including required/optional fields that are missing.

$ from HEK import HER_Event

$ event = HER_Event("FL")
$ event.required["OBS_OBSERVATORY"] = 'TRACE'
$ event.required["OBS_INSTRUMENT"] = 'TRACE'
$ event.required["OBS_CHANNELID"] = 'TRACE 171'
$ event.required["OBS_MEANWAVEL"] = '171'
$ event.required["OBS_WAVELUNIT"] = 'Angstroms'
$ event.required["FRM_NAME"] = 'Karel Schrijver'
$ event.required["FRM_IDENTIFIER"] = 'Karel Schrijver'
$ event.required["FRM_INSTITUTE"] ='LMSAL'
$ event.required["FRM_HUMANFLAG"] = 'yes'
$ event.required["FRM_PARAMSET"] = 'n/a'
$ event.required["FRM_DATERUN"] = '2007/01/03 12:00:00'
$ event.required["FRM_CONTACT"] = 'schryver at lmsal dot com'
$ event.required["EVENT_STARTTIME"] = '2006/10/10 23:45:13'
$ event.required["EVENT_PEAKTIME"] = '2006/10/10 23:47:54'
$ event.required["EVENT_ENDTIME"] = '2006/10/10 23:55:20'
$ event.required["EVENT_COORDSYS"] = 'UTC-HPC-TOPO'
$ event.required["EVENT_COORDUNIT"] = 'arcsec'
$ event.required["EVENT_COORD1"] = '-400'
$ event.required["EVENT_COORD2"] = '300'
$ event.required["EVENT_C1ERROR"] = '4'
$ event.required["EVENT_C2ERROR"] = '4'
$ event.required["BOUNDBOX_C1LL"] = '-440'
$ event.required["BOUNDBOX_C2LL"] = '260'
$ event.required["BOUNDBOX_C1UR"] = '-360'
$ event.required["BOUNDBOX_C2UR"] = '340'
$ event.Description = 'My first flare'

$ print(event)

$ event.exportEvent()

>> Event codes <<
AR - Active Region
CE - CME
CD - Coronal Dimming
CW - Coronal Wave
FI - Filament
FE - Filament Eruption
FA - Filament Activation
FL - Flare
LP - Loop
OS - Oscillation
SS - Sunspot
EF - Emerging Flux
CJ - Coronal Jet
PG - Plage
OT - Other
NR - Nothing Reported
SG - Sigmoid
SP - Spray Surge
CR - Coronal Rain
CC - Coronal Cavity
ER - Eruption

Note: The "citations" section has not yet been implemented. All references are instead placed at the bottom of the xml tree root.
