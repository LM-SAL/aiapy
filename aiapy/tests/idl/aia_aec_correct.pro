pro AIA_AEC_CORRECT, hdrs, $
  dnbias = dnbias, exptimebias = exptimebias, shortonly = shortonly, $
  satscale = satscale, satshift = satshift, $
  win = win, pstop = pstop, $
  fixcounts_on = fixcounts_on, fixcounts_off = fixcounts_off, noplot = noplot
  compile_opt idl2

  ; aiapy: added the fixcounts_on/fixcounts_off outputs and /noplot for the IDL comparison test.

  ;+
  ;
  ; Makes a plot showing the correction to the count rate for short exposures
  ;
  ; INPUT:
  ;	hdrs	-	A structure array of image headers
  ;
  ; KEYWORDS:
  ;	dnbias	-	Number of DN in dark offset
  ;	exptimebias -	Error in exposure time (in msec)
  ;	/shortonly	-	If set, then the exposure time error is only applied to short
  ;					(narrow-slit) exposures. This should be the default but isn't
  ;	satscale	-	Scaling between NSATPIX and DN (DN of ghosting per NSATPIX)
  ;					Try something like 18000
  ;	satshift	-	Time lag between saturation and ghosting.
  ;
  ;-

  ; Here's a sample setup script:
  ; SSW_JSOC_TIME2DATA, '2014-09-08T23:00', '2014-09-09T04:00', drms, files, /files, ds = 'aia.lev1', wave = 131
  ; SSW_JSOC_TIME2DATA, '2014-09-10T17:00', '2014-09-10T21:00', drms, files, /files, ds = 'aia.lev1_nrt2', wave = 131
  ; READ_SDO, files, hdrs, data
  ; data_corner = data[0:63, 0:63, *]
  ; cmean = TOTAL(TOTAL(data_corner, 1), 1) / 64. / 64.
  ; PLOT, cmean
  ; dnbias = MEAN(cmean)
  ; goodpix = WHERE(drms.exptime ne 0 and drms.exptime / drms.aimgshce * 1000 gt 0.8)
  ; AIA_AEC_CORRECT, drms[goodpix], dnbias = 0-dnbias, /short, exptime = 0.0032, satscale = 18e3

  if n_elements(dnbias) eq 0 then dnbias = 0.0
  dnbias_str = string('X = ', dnbias, form = '(a4,f7.4)')
  if n_elements(exptimebias) eq 0 then exptimebias = 0.0
  exptimebias_str = string('Y = ', exptimebias, form = '(a4,f7.4)')
  if n_elements(satscale) eq 0 then satscale = 0
  biasstr = [dnbias_str, exptimebias_str]
  shortonly = keyword_set(shortonly)

  i_aecon = where(hdrs.aectype gt 1, numon)
  t_on = hdrs[i_aecon].t_obs
  tai_on = ANYTIM2TAI(t_on)
  mean_on = hdrs[i_aecon].datamean
  exptime_on = hdrs[i_aecon].exptime

  i_aecoff = where(hdrs.aectype le 1, numoff)
  t_off = hdrs[i_aecoff].t_obs
  tai_off = ANYTIM2TAI(t_off)
  mean_off = hdrs[i_aecoff].datamean
  exptime_off = hdrs[i_aecoff].exptime

  ; If the /shortonly keyword is set, then only apply the exposure time correction
  ; to short (narrow slit) exposures
  if shortonly eq 1 then begin
    shorts = where(exptime_on lt 75)
    exptimebias_arr = exptime_on * 0.
    exptimebias_arr[shorts] = exptimebias
    corr_exptime_on = (exptime_on - exptimebias_arr)
    corr_exptime_off = exptime_off
  endif else begin
    corr_exptime_on = (exptime_on - exptimebias)
    corr_exptime_off = (exptime_off - exptimebias)
  endelse

  rawcounts_on = mean_on / exptime_on
  rawcounts_off = mean_off / exptime_off
  fixcounts_on = (mean_on - dnbias) / corr_exptime_on
  fixcounts_off = (mean_off - dnbias) / corr_exptime_off

  ; If the satscale is set to something, then apply a saturation correction
  if satscale ne 0 then begin
    if n_elements(satshift) eq 0 then satshift = -12
    satcounts = hdrs[i_aecoff].nsatpix / float(satscale)
    satdelta_on = interpol(satcounts, tai_off + satshift, tai_on)
    fixcounts_on = fixcounts_on - satdelta_on
    satcounts_off = hdrs[i_aecon].nsatpix / float(satscale)
    satdelta_off = interpol(satcounts_off, tai_on + satshift, tai_off)
    fixcounts_off = fixcounts_off - satdelta_off
    satscale_str = string('Satscale = ', satscale, form = '(a11, i6)')
    satshift_str = string('Satshift = ', satshift, form = '(a11, f6.1)')
    biasstr = [biasstr, satscale_str, satshift_str]
    ytitle = '(Datamean - X)/(Exptime - Y) - NSATPIX_shift [DN/s/pix]'
  endif else begin
    ytitle = '(Datamean - X)/(Exptime - Y) [DN/s/pix]'
  endelse
  if keyword_set(noplot) then RETURN
  ifixcounts_off = interpol(fixcounts_off, tai_off, tai_on)
  irawcounts_off = interpol(rawcounts_off, tai_off, tai_on)
  flaretime = where(smooth(abs(rawcounts_on - irawcounts_off), 10) gt 0.01, numflare)
  rmsraw = mean(abs(rawcounts_on[flaretime] - irawcounts_off[flaretime]))
  rmsfix = mean(abs(fixcounts_on[flaretime] - ifixcounts_off[flaretime]))
  biasstr = [biasstr, 'Shortonly = ' + strtrim(shortonly, 2)]

  if n_elements(win) gt 0 then WDEF, win, 1000, 800
  UTPLOT, t_on, fixcounts_on, psym = -4, /ynoz, chars = 1.5, /xstyle, $
    ytitle = ytitle, title = 'Full disk signal', $
    pos = [0.07, 0.3, 0.95, 0.95]
  OUTPLOT, t_off, fixcounts_off, psym = -4, col = 2
  OUTPLOT, t_on, rawcounts_on, psym = -4, col = 6
  LEGEND, pos = 12, chars = 2, col = [255, 2, 6], ['AEC off', 'AEC on', 'AEC on (raw)']
  LEGEND, pos = 10, chars = 1.5, biasstr, /nolines
  UTPLOT, t_on, rawcounts_on - irawcounts_off, /noerase, chars = 1.5, /xstyle, $
    ytit = 'AEC On - AEC off', $
    pos = [0.07, 0.01, 0.95, 0.29]
  OUTPLOT, t_on, fixcounts_on - ifixcounts_off, col = 2
  LEGEND, pos = 12, chars = 2, col = [255, 2], [string('Raw RMS = ', rmsraw, form = '(a10,f6.4)'), $
    string('Corr. RMS = ', rmsfix, form = '(a12,f6.4)')]
  PB_WIN2PNG, !desk + 'scratch/aia_aec_correct.png'

  if keyword_set(pstop) then stop
end
