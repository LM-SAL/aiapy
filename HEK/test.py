from HEK import HER_Event

event = HER_Event("FL")
event.required["OBS_OBSERVATORY"] = 'TRACE'
event.required["OBS_INSTRUMENT"] = 'TRACE'
event.required["OBS_CHANNELID"] = 'TRACE 171'
event.required["OBS_MEANWAVEL"] = '171'
event.required["OBS_WAVELUNIT"] = 'Angstroms'
event.required["FRM_NAME"] = 'Karel Schrijver'
event.required["FRM_IDENTIFIER"] = 'Karel Schrijver'
event.required["FRM_INSTITUTE"] ='LMSAL'
event.required["FRM_HUMANFLAG"] = 'yes'
event.required["FRM_PARAMSET"] = 'n/a'
event.required["FRM_DATERUN"] = '2007/01/03 12:00:00'
event.required["FRM_CONTACT"] = 'schryver at lmsal dot com'
event.required["EVENT_STARTTIME"] = '2006/10/10 23:45:13'
event.required["EVENT_PEAKTIME"] = '2006/10/10 23:47:54'
event.required["EVENT_ENDTIME"] = '2006/10/10 23:55:20'
event.required["EVENT_COORDSYS"] = 'UTC-HPC-TOPO'
event.required["EVENT_COORDUNIT"] = 'arcsec'
event.required["EVENT_COORD1"] = '-400'
event.required["EVENT_COORD2"] = '300'
event.required["EVENT_C1ERROR"] = '4'
event.required["EVENT_C2ERROR"] = '4'
event.required["BOUNDBOX_C1LL"] = '-440'
event.required["BOUNDBOX_C2LL"] = '260'
event.required["BOUNDBOX_C1UR"] = '-360'
event.required["BOUNDBOX_C2UR"] = '340'
event.Description="My first flare"

print (event)

event.exportEvent()
