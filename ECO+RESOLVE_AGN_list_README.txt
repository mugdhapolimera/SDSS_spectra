ECO+RESOLVE_AGN_list README
July 7th 2021
This file is a master list of ECO and RESOLVE AGN identified till date. 

Column name 	Units		Description
Name				ECO/RESOLVE identifier. Caution: Some RESOLVE galaxies 
				may also be in ECO, but are identified only by their 
				RESOLVE names
radeg		deg		RA in degrees
dedeg		deg		Dec in degrees
vhel		km/s		Heliocentric recession velocity
logmstar			log(stellar mass/M_sun)
logmgas				log(gas mass/M_sun)
logmh				log(halo mass/M_sun) from halo abundance matching (as in 
				RESOLVE/ECO public catalogs)
gmag		mag		Apparent SDSS g-band magnitude (without 
				foreground extinction corrections)
rmag		mag		Apparent SDSS r-band magnitude (without 
				foreground extinction corrections)
oiii_5007_flux	10^-17 ergs/s/cm^2	[O III] 5007A flux from MPA-JHU SDSS catalog 
				with extinction correction
oiii_5007_flux_err	10^-17 ergs/s/cm^2	Error on [O III] 5007A flux from MPA-JHU SDSS catalog
o3lum		ergs/s		Luminosity from oiii_5007_flux 
nuvmag				Apparent GALEX NUV-band magnitude (without 
				foreground extinction corrections)
rband_censb	mag/arcsec^2	Central surface brightness (within 3" from center) in 
				r-band 
obstel				Intended/used telescope for observations. 
				N/A - no observations intended; SALT - SALT long-slit 
				emission line setup; SALTblue - SALT long-slit abs. line
				setup; SOAREM - SOAR slicer emission line setup; SOARABS - 
				SOAR slicer abs. line setup; GEMEM - Gemini IFU emission
				line setup
obsstat				Observation status in RESOLVE; notobs - not observed; 
				bluedone/ - obs done in abs. line setup; reddone/ - obs 
				done in emission line setup; broaddone/ - obs done in
				SOAR broad setup. More than one setup can be done (e.g.,
				bluedone/broaddone)
broadline			True flag if match available in Liu et al. 2019 broadline
				AGN catalog
agntype				Type of AGN identified. sfingagn - SF-AGN from Polimera 
				et al. (in prep.), SF in BPT (NII plot) but AGN in VO87
				(SII and OI plots); bptagn - Seyfert or LINER in BPT;
				bptcomposite - Composite in BPT; xrayagn - AGN from 
				L_xray-SFR relationship by Ranalli+ 03; midiragn - AGN by
				any mid-IR colour criteria (Jarrett+12, Satyapal+14, Stern+11), 
				midiragn only available for RESOLVE as of now. 
ifusource			'manga'/'sami' in case of crossmatches, 'none' otherwise
