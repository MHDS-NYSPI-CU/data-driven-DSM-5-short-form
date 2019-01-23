%macro DSF	(	file=, 				/*REQUIRED - file reference for where data is located*/
							dsn=,				/* REQUIRED - person-level dataset provided by investigator with DSM criteria, chosen demographics, weighting */
							criteriavars=, 		/* REQUIRED - list of criteria variables in datasets */
							dxvar=,				/* REQUIRED - dsm diagnosis variable */
							wtvar=, 				/* OPTIONAL - weight var if used in dataset - weight will be used to obtain estimates */
							DataPrep=No,	/* OPTIONAL - if No then Data section does NOT run, must turn on --> Yes (DEFAULT=No) */
							Corr=No,			/* OPTIONAL - if No then Correlation estimates NOT created, must turn on --> Yes (DEFAULT=No) */
							Prev=No,			/* OPTIONAL - if No then Prevalence estimates NOT created, must turn on --> Yes (DEFAULT=No) */
							SeSp=No,			/* OPTIONAL - if No then Sens/Spec estimates NOT created, must turn on --> Yes (DEFAULT=No) */
							bsamp=50, 		/* OPTIONAL - number of bootstrap samples to be created  (DEFAULT=50) */
							pdiff=0.05, 		/*OPTIONAL - DSFs within % range of current rule (DEFAULT=5%=0.05) */
							diagacc=0.90, 	/* OPTIONAL - diagnostic accuracy threshold for Sensitivity/Specificity estimates >= (DEFAULT=0.90) */
							DIFF=No,   		/* OPTIONAL - if No then Differential Functioning by Covariates NOT assessed, must turn on --> Yes (DEFAULT=No) */
							covars=				/* REQUIRED - list of covariates to assess Differential Functioning by Covariates */
								
);

options ls=80 ps=55 nodate nofmterr mprint mlogic symbolgen MergeNoBy=error MSGLevel = I compress=yes;
libname dsf "&file";

%if "&DataPrep"="Yes" %then %do;

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Starting Data Prep Section - Creates Data for Algorithm.;
%put;
options nonotes;

********************************************************************************************************************************************************************;
*	(A) DATA PREP																																													;
********************************************************************************************************************************************************************;

data temp_FinalCombo1;run;
data temp_FinalCombo1;
	length CriteriaList $200. CriteriaNumList $20.;
	set temp_FinalCombo1;
	p=0;
	CriteriaList="test";
	CriteriaNumList="x";
run;

*******************************************************;
*	CREATE BOOTSTRAP SAMPLES				;
*******************************************************;

* Create a sequential patid, numbered 1 to n using investigator provided dataset;
data temp_BootIDData0;
	set dsf.&dsn;
	patid + 1;
run;

* Obtain n in sample ds;
proc sql;
	select count(patid) into :n trimmed
	from temp_BootIDData0;
quit;

* Generate # boostrap samples of random nesarc ids;
*		Dataset BootSampID provides IDs to link original data with bootstrap data									;
*		BootSamp = 1 to # of bootstrap samples you want (with n cases in each sample)						;
*		BootSampPatid = 1 to n sequential number within each BootSamp											;
*		RandPatid = creates a random patid for each BootSampPatid (this will be linked to original data)	;
data BootSampID;
	do BootSamp=1 to &bsamp;
		do BootSampPatid=1 to &n;
			RandPatid=ceil(&n*ranuni(1));
			output;
		end;
	end;
run;

*Link Random PATID (RandPatid) from Bootstrap sample (BootSampID) to the original data	;
proc sql;
	create table temp_BootIDData1 as
	select distinct a.*, b.BootSamp, b.BootSampPatid
	from temp_BootIDData0 a inner join BootSampID b
		on a.patid=b.RandPatid
		order by BootSamp, BootSampPatid;
quit;

*******************************************************;
*	CREATE MACRO CALLS							;
*******************************************************;

* Create macro call for criteria var names in with parenthesis " ";
*	&criteriaparen= "criteria1" "criteria2"... "criteriaN"					;
data temp_criteria;run;
data temp_criteria;
	length tempcorrs1 tempcorrs2 $500.;
	set temp_criteria;

	tempcorrs1=upcase("&criteriavars");
	tempcorrs2=left(trim('"'||strip(TRANWRD(trim(tempcorrs1),trim(' '),'" "'))||'"'));
	call symput("criteriaparen",strip(tempcorrs2));
run;

* Create macro calls for criteria vars ;
*	&Crt1=criteria variable name (varname for criteria 1)	;
*	&Crtn1=number for criteria (1)									;
*	&nCrt = total number of criteria								;
proc contents data=dsf.&dsn memtype=data 
	out=CriteriaVars0 (keep=NAME LABEL VARNUM where=(upcase(NAME) in (&criteriaparen))) 
	noprint;
run;
data CriteriaVars;
	set CriteriaVars0;
	VARNUMx = input(compress(NAME, "","A"),best.);
run;
proc sort data=CriteriaVars; by VARNUMx;run;
data _null_;
	set CriteriaVars;
	call symputx(compress("Crt"||put(_n_,best.)), NAME);
	call symputx (compress("Crtn"||put(_n_,best.)),put(_n_,best.));
run;
%global nCrt;
proc sql noprint;
	select count(NAME) into :nCrt trimmed
	from CriteriaVars;
quit;

*******************************************************;
*	CREATE CRITERIA COMBINATIONS			;
*******************************************************;

*macro creates all possible combinations of criteria (subscales)	;
*creates Final_Varlist_0 - all criteria combos								;
%macro mix(Combo=,j=1,start=1,xmod=, zmod=);
	%do i&j=&start %to &nCrt;
		%let model=&xmod &&&&Crt&&i&j;
		%let numlist=&zmod &&&&Crtn&&i&j;
		%**Discriminat analysis***;
			data temp_FinalCombo2;
				length CriteriaList $200. CriteriaNumList $20.;
				set temp_FinalCombo1;
				CriteriaList = "&model";
				CriteriaNumList=compress("&numlist");
				p = &j; *p = length of criteria used (total number of variables used);
			run;
		%**Append information to ouput dataset***;
			data FinalCombo;
				set 	FinalCombo
						temp_FinalCombo2;
				if CriteriaList = "" then delete;
			run;
		%if &j<&Combo %then 
			%mix(Combo=&Combo,j=%eval(&j+1), start=%eval(&&i&j+1),xmod=&model, zmod=&numlist);
		%end;
%mend mix;
*initiate datasets;
data FinalCombo;run;
%mix(Combo=&nCrt); *all possible combinations;

*Turn Log info back on;
options notes mprint source source2 symbolgen mlogic;

*Add in binary for criteria;
data FinalCombox;
	set FinalCombo;
	%macro criteria;
		%do i=1 %to &nCrt;
			if index(CriteriaList,"&&Crt&i ")>0 then C&i=1;
			else C&i=0;
		%end;
	%mend criteria;
	%criteria;
run;

*******************************************************;
*	CREATE COMBO MACRO CALLS				;
*******************************************************;

*Create macro calls for all the possible combinations of criteria 	;
*	&nCombo = total number of possible combinations of criteria 	;
*	&Combo# = list of criteria variables for specified combination	;

*Make GLOBAL Macros;
%global nCombo;
data _null_;
	set FinalCombox;
	call symputx (compress("Combo"||put(_n_,best.)),CriteriaList);
run;
proc sql noprint;
	select count(CriteriaList) into :nCombo trimmed
	from FinalCombox;
quit;
data FinalCombox;
	set FinalCombox;
	CRule=_n_;
run;

*******************************************************;
*	CREATE SUBSCALES								;
*******************************************************;

*Create Subscales by summing all possible combinations of criteria;
data BootIDData1;
	set temp_BootIDData1;
	%macro sum;
		%do i=1 %to &nCombo;
			Subscale&i = sum(of &&Combo&i);
		%end;
	%mend;
	%sum;
run;

*DELETE TEMP DS;
proc datasets lib = work;
	delete temp:;
quit;
run;

*OUTPUT PERM DATASETS;
data dsf.BootIDData1; set BootIDData1;run;
data dsf.BootSampID; set BootSampID;run;
data dsf.CriteriaVars; set CriteriaVars;run;
data dsf.FinalCombo; set FinalCombox;run;

%end; *End DataPrep section;

********************************************************************************************************************************************************************;
*	(B) CREATE ESTIMATES AND 95% CIs																																				;
********************************************************************************************************************************************************************;

*******************************************************;
* B1: CORRELATION W/ TOTAL						;
*******************************************************;

%if "&Corr"="Yes" %then %do;

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Correlation with Total - Creates subscales for all possible combinations and correlates with total.;
%put;
options nonotes;

*Call in Previously Created Data;
data CriteriaVars; set dsf.CriteriaVars;run;
data FinalCombo; set dsf.FinalCombo;run;

%global nCrt nCombo;
proc sql noprint;
	select count(NAME) into :nCrt trimmed
	from CriteriaVars;
quit;
proc sql noprint;
	select count(CriteriaList) into :nCombo trimmed
	from FinalCombo;
quit;

*Estimate corrleation with total;
ods select none;
ods output PearsonCorr=temp_CorrData0;
proc corr data=dsf.BootIDData1 pearson;
	%if %length(&wtvar)>0 %then %do;
		weight &wtvar;
	%end;
	by BootSamp;
	var Subscale&nCrt; *subscale using all criteria - this is TOTAL;
	with Subscale1-Subscale&nCombo; *go through all possible subscales;
run;
ods select all;

*Finalize Boot-Subscale level data;
data temp_CorrData1 (keep=BootSamp PCorr&nCrt CRule);
	set temp_CorrData0 (rename=(Subscale&nCrt=PCorr&nCrt));

	*Pull Subscale info;
	%macro length1;
		%do i=9 %to 12;
			if length(Variable)=&i then CRule=input(substr(Variable,9,%eval(&i-8)),best.);
		%end;
	%mend length1;
	%length1;
run;

*Estimate 95% CIs;
proc sort tagsort data=temp_CorrData1; by CRule;run;
proc means data=temp_CorrData1 noprint;
	by CRule;
	vars PCorr&nCrt;
	output out=temp_CorrCI0 mean=PCorr&nCrt.Mean std=PCorr&nCrt.Std;
run;
data temp_CorrCI1;
	set temp_CorrCI0 (drop=_:);

	*95% CI;
	PCorr&nCrt.CI95L = PCorr&nCrt.Mean - (1.96*PCorr&nCrt.Std);
	PCorr&nCrt.CI95U = PCorr&nCrt.Mean + (1.96*PCorr&nCrt.Std);
run;

*Add Correlation Estimates & 95% CIs to final Subscale Level dataset;
data 	CorrCIData;
	merge	FinalCombo (in=a)
				temp_CorrCI1 (drop=PCorr&nCrt.Std);
	by CRule;
run;

*Identify Top Performing Subscales (within p);
data FinalCorrCI;run;

%do i=1 %to %eval(&nCrt-1); *runs through all sizes of combinations (p-1) excludes full version;

	proc sort data=CorrCIData out=temp_P&i._0 (where=(p=&i)); by descending PCorr&nCrt.Mean;run;
	*pull top performer - create macro for 95% LB;
	data _null_;
		set temp_P&i._0 (firstobs=1 obs=1);
		call symputx("cutpt",PCorr&nCrt.CI95L);
	run;
	%put &cutpt;
	*pull any p short-forms >= 95% LB from top performer;
	data temp_P&i;
		set temp_P&i._0 (where = (PCorr&nCrt.CI95U >= &cutpt));
	run;
	*add top short-forms to final dataset with all p;
	data FinalCorrCI;
		set 	FinalCorrCI
				temp_P&i;
		if p = . then delete;
	run;
	*delete ds;
	proc datasets lib = work;
		delete temp_P&i.:;
	quit;
	run;
%end;

*DELETE TEMP DS;
proc datasets lib = work;
	delete temp:;
quit;
run;

*******************************************************;
* CREATE DSFs FROM TOP CORRS				;
*******************************************************;

*Create Macro Calls to create DSFs (used in Prev too)	;
*	cCrt1 = varname of criteria included in combo 1		;
*	cP1 = size of combo 1											;
*	cDSF = name for DSF -- DSF[2047] 						;
*	cNum = 2047 reference number								;
*	cN = total number of subscales selected 					;
%global cN;
data FinalCorrCI;
	set FinalCorrCI;
	CorrNum=_n_;
	call symputx (compress("cCrt"||put(CorrNum,best.)),CriteriaList);
	call symputx (compress("cP"||put(CorrNum,best.)),put(p,best.));
	call symputx (compress("cDSF"||put(CorrNum,best.)),compress("DSF"||put(CRule,best.)));
	call symputx (compress("cNum"||put(CorrNum,best.)),CRule);
run;
proc sql noprint;
	select count(CriteriaList) into :cN trimmed
	from FinalCorrCI;
quit;

*Sum all combinations at the ID-Boot Level & create dichotomous rules;
data BootIDData2;
	set dsf.BootIDData1 (drop=Subscale:);
	%macro sum;
		%do i=1 %to &cN;
			&&cDSF&i = sum(of &&cCrt&i);

			*Create dichotomous rules;
			%do m=1 %to &&cP&i;
				if &&cDSF&i >= &m then &&cDSF&i.._m&m = 1;
				else if 0 <= &&cDSF&i < &m then &&cDSF&i.._m&m = 0;
			%end;					
		%end;
	%mend sum;
	%sum;
run;
proc sort tagsort data=BootIDData2; by BootSamp;run;

*DELETE TEMP DS;
proc datasets lib = work;
	delete temp: BootIDData1;
quit;
run;

*OUTPUT PERM DATASETS;
data dsf.BootIDData2; set BootIDData2;run;
data dsf.CorrCIData; set CorrCIData;run;
data dsf.FinalCorrCI; set FinalCorrCI;run;

%end; *End Corr section;

*******************************************************;
* B2: COMPARE PREV OF DSFs VS CURRENT RULE	;
*******************************************************;

%if "&Prev"="Yes" %then %do;

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Prevalence - Create DSFs and compare prevalence to prevalence using current rule.;
%put;
options nonotes;

*Create Macro Call for AUD prevalence in sample;
*	DSMp = prevalence of DSM disorder in sample;
proc means data=dsf.BootIDData2;
	by BootSamp;
	%if %length(&wtvar)>0 %then %do;
		weight &wtvar;
	%end;
	vars &dxvar;
	output out=tempDSMPrev0 mean=prev;
run;
proc means data=tempDSMPrev0;
	vars prev;
	output out=tempDSMPrev1 mean=mDSMPrev std=stdDSMPrev;
run;
data DSMPrev;
	set tempDSMPrev1 (drop=_:);
	mDSMPrevCI95LB = mDSMPrev - (1.96*stdDSMPrev);
	mDSMPrevCI95UB = mDSMPrev + (1.96*stdDSMPrev);
	m=2;
	p=11;
run;
proc sql noprint;
	select mDSMPrev, mDSMPrevCI95LB, mDSMPrevCI95UB into :dsmPrev trimmed, :dsmPrev95L trimmed, :dsmPrev95U trimmed
	from DSMPrev;
quit;

**************************************;
*Create Prevalence Estimates & Differences for All DSFs (from Corr);
**************************************;

*Re-use macros;
data FinalCorrCI;
	set dsf.FinalCorrCI;
	CorrNum=_n_;
	call symputx (compress("cCrt"||put(CorrNum,best.)),CriteriaList);
	call symputx (compress("cP"||put(CorrNum,best.)),put(p,best.));
	call symputx (compress("cDSF"||put(CorrNum,best.)),compress("DSF"||put(CRule,best.)));
	call symputx (compress("cNum"||put(CorrNum,best.)),CRule);
run;
proc sql noprint;
	select count(CriteriaList) into :cN trimmed
	from FinalCorrCI;
quit;

*initiate ds;
data temp_xDiff;run;
data temp_xPrev;run;

*Go through each DSF created from the subscales selected from CORR Section;
%do i=1 %to &cN; 

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Calculating Prevalence Estimates (&i of &cN at the subscale level - broken down further within macro);
%put;
options nonotes;

	*Create m for each subscale (m=1 to p);
	%do m=1 %to &&cP&i; *the p for the specified subscale;

		**************************************;
		*Difference (DSF prev - AUD prev) + 95 % CIS;
		**************************************;
		proc means data=dsf.BootIDData2;
			by BootSamp;
			%if %length(&wtvar)>0 %then %do;
				weight &wtvar;
			%end;
			vars &&cDSF&i.._m&m;
			output out=temp_Prev0 mean=prev;
		run;
		data temp_Diff0;
			set temp_Prev0;
			DSM=&dsmPrev;
			Diff=(prev-DSM);
		run;
		*obtain 95% CIs for Diff in all 50 bootstrap samples;
		proc means data=temp_Diff0;
			vars Diff;
			output out=temp_Diff1 mean=mPDiff std=stdPDiff;
		run;
		data temp_Diff2;
			length dsf $10.;
			set temp_Diff1 (drop=_:);
			mPDiff95L = mPDiff - (1.96*stdPDiff);
			mPDiff95U = mPDiff + (1.96*stdPDiff);
			dsf=left(trim(compress("&&cDSF&i.._m&m")));
			m=&m;
			p=&&cP&i;
			CRule=&&cNum&i;
		run;
		data temp_xDiff;
			set 	temp_xDiff
					temp_Diff2;
			if dsf = "" then delete;
		run;
		*obtain 95% CIs for Prev in all 50 bootstrap samples;
		proc means data=temp_Prev0;
			vars prev;
			output out=temp_Prev1 mean=mPrev std=stdPrev;
		run;
		data temp_Prev2;
			length dsf $10.;
			set temp_Prev1 (drop=_:);
			mPrev95L = mPrev - (1.96*stdPrev);
			mPrev95U = mPrev + (1.96*stdPrev);
			dsf=left(trim(compress("&&cDSF&i.._m&m")));
		run;
		data temp_xPrev;
			set 	temp_xPrev
					temp_Prev2;
			if dsf = "" then delete;
		run;
	%end; *end m;
%end; *end through all dsf;

*Turn Log info back on;
options notes mprint source source2 symbolgen mlogic;

*Final Dataset with Diff + 95% CIs AND Prev + 95% CIs;
proc sort data=temp_xPrev;by dsf;run;
proc sort data=temp_xDiff;by dsf;run;
data temp_FinalPrevCI0;
	merge	temp_xPrev (in=a)
				temp_xDiff (in=b);
	by dsf;
	if a or b;
	mPDiffAbs=abs(mPDiff);
run;

*add in rule info & sort by abs diff to get correct order (we do not plot this);
proc sort data=temp_FinalPrevCI0; by CRule p;run;
proc sort data=dsf.FinalCombo out=FinalCombo; by CRule p;run;
data temp_FinalPrevCI1;
	merge	temp_FinalPrevCI0 (in=a)
				FinalCombo;
	by CRule p;
	if a;
run;
proc sort data=temp_FinalPrevCI1; by mPDiffAbs;run;
data 	FinalPrevCI
		PrevCIData;
	set temp_FinalPrevCI1;
	by mPDiffAbs;
	*Within % of DSM Prev;
	if - &pdiff <= mPDiff95L <= &pdiff then PrevKeep=1;
	if - &pdiff  <= mPDiff95U <= &pdiff then PrevKeep=1;
	output PrevCIData;
	if PrevKeep=1 then output FinalPrevCI;
run;

*DELETE TEMP DS;
proc datasets lib = work;
	delete temp:;
quit;
run;

*OUTPUT PERM DATASETS;
data dsf.PrevCIData; set PrevCIData; run;
data dsf.FinalCorrCI; set FinalCorrCI;run;
data dsf.FinalPrevCI; set FinalPrevCI;run;
data dsf.DSMPrev; set DSMPrev;run;

%end; *End Prev Section;

%if "&SeSp"="Yes" %then %do;

*******************************************************;
* B3: SENSITIVITY & SPECIFICITY 				;
*******************************************************;

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Sensitivity/Specificity - Assesses diagnostic accuracy for each DSF compared to current rule.;
%put;
options nonotes;

*Create Macro Calls to create DSFs (from Prev)			;
*	pP1 = size of dsf 1													;
*	pDSF = name for DSF -- DSF[2047]_m#					;
*	pNum = 2047 reference number								;
*	pN = total number of subscales selected 					;

data FinalPrevCI;
	set dsf.FinalPrevCI;
	PrevNum=_n_;
	call symputx (compress("pP"||put(PrevNum,best.)),put(p,best.));
	call symputx (compress("pM"||put(PrevNum,best.)),put(m,best.));
	call symputx (compress("pDSF"||put(PrevNum,best.)),dsf);
	call symputx (compress("pNum"||put(PrevNum,best.)),CRule);
run;
proc sql noprint;
	select count(dsf) into :pN trimmed
	from FinalPrevCI;
quit;

**************************************;
*Go through all the optimal DSFs to create SENS/SPEC;
**************************************;
data BootIDData2;
	set dsf.BootIDData2;
run;

*initiate final Sens/Spec dataset;
data SSCIData0;run;

%do j=1 %to &pN;

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Calculating Sensitivity/Specificity Estimates (&j of &pN);
%put;
options nonotes;

	ods html close;
	ods preferences;
	ods html newfile=none;
	*Calculate Sensitivity & Specificity - all boots;
	proc freq data=BootIDData2;
		by BootSamp;
		weight &wtvar;
		tables &&pDSF&j * &dxvar;
		ods output crosstabfreqs=temp_Freq0 (keep = BootSamp Table &&pDSF&j &dxvar Frequency ColPercent);
	run;

	*pull of sensitivity and specificity;
	data temp_Freq1_sens (keep=BootSamp Sensitivity TP)
			temp_Freq1_spec (keep = BootSamp Specificity: TN)
			temp_Freq1_FP (keep = BootSamp FP)
			temp_Freq1_FN (keep = BootSamp FN);
		set temp_Freq0;
		*Sensitivity + TP;
		if &&pDSF&j = 1 and &dxvar = 1 then do;
			Sensitivity = round((ColPercent/100),.00001);
			TP=Frequency;
			output temp_Freq1_sens;
		end;
		*Specificity + TN;
		if &&pDSF&j = 0 and &dxvar = 0 then do;
			Specificity = round((ColPercent/100),.00001);
			Specificity1Minus = 1-Specificity;
			TN=Frequency;
			output temp_Freq1_spec;
		end;
		*FP;
		if &&pDSF&j = 1 and &dxvar = 0 then do;
			FP = Frequency;
			output temp_Freq1_FP;
		end;
		*FN;
		if &&pDSF&j = 0 and &dxvar = 1 then do;
			FN = Frequency;
			output temp_Freq1_FN;
		end;
	run;
	*merge all stats together for 1 row per bootstrap (should be 50 total);
	data temp_Freq2;
		merge	temp_Freq1_sens (in=a)
					temp_Freq1_spec
					temp_Freq1_FP
					temp_Freq1_FN;
		by BootSamp;
		if a;
	run;

	*mean of all sens/spec over all boots;
	%macro sens_spec (type=, short=);
		proc means data=temp_Freq2;
			vars &type;
			output out=temp_Freq2_&short.0 mean=m&short std=std&short;
		run;
		*Creat 95% CIs;
		data temp_Freq2_&short;
			length DSF $10.;
			set temp_Freq2_&short.0 (drop=_:);
			m&short.95LB = m&short - (1.96*std&short);
			m&short.95UB = m&short + (1.96*std&short);
			DSF=trim(left(compress("&&pDSF&j")));
		run;
	%mend sens_spec;
	%sens_spec (type=Sensitivity, short=Sens);
	%sens_spec (type=Specificity, short=Spec);

	*Merge final meanSens and meanSpec into final dataset (DSF level);
	data temp_SSCIData;
		merge	temp_Freq2_sens
					temp_Freq2_spec;
		by DSF;
		CRule=&&pNum&j;
		p=&&pP&j;
		m=&&pM&j;
	run;
	*add to final SensSpec dataset;
	data SSCIData0;
		set SSCIData0
			temp_SSCIData;
		if DSF=" " then delete;
	run;
	*delete freq datasets;
	proc datasets lib = work;
		delete temp_:;
	quit;
	run;
%end;

*Turn Log info back on;
options notes mprint source source2 symbolgen mlogic;

*Add in dsf info;
proc sort data=SSCIData0; by CRule p;run;
proc sort data=dsf.FinalCombo out=FinalCombo; by CRule p;run;
data SSCIData
		FinalSSCI;
	merge	SSCIData0 (in=a)
				FinalCombo (in=b);
	by CRule p;
	if a;

	*temp added - remove once i add 1-spec to macro;
	mSpec1=1-mSpec;

	*Significantly >= 0.90 for Sens & Spec;
	if mSens >=&diagacc or mSens95LB >= &diagacc or mSens95UB >= &diagacc then KeepSens=1;
	if mSpec >=&diagacc or mSpec95LB >= &diagacc or mSpec95UB >= &diagacc then KeepSpec=1;

	output SSCIData;
	if KeepSens=1 and KeepSpec=1 then output FinalSSCI;
run;

*DELETE TEMP DS;
proc datasets lib = work;
	delete SSCIData0;
quit;
run;

*OUTPUT PERM DATASETS;
data dsf.SSCIData;set SSCIData;run;
data dsf.FinalSSCI;
	set FinalSSCI;
	SSNum=_n_;
run;

%end; *End Sens/Spec Section;

%if "&DIFF"="Yes" %then %do;

*******************************************************;
* B4: DIFFERENTIAL FUNCTIONING BY COVARIATES				;
*******************************************************;

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Differential Functional by Covariates - Assess DIFF for selected covaraites.;
%put;
options nonotes;

*Create Macro Calls to go through DSFs from SS			;
*	pP1 = size of dsf 1													;
*	pDSF = name for DSF -- DSF[2047]_m#					;
*	pNum = 2047 reference number								;
*	pN = total number of subscales selected 					;
data FinalSSCI;
	set dsf.FinalSSCI;
	SSNum=_n_;
	call symputx (compress("sP"||put(SSNum,best.)),put(p,best.));
	call symputx (compress("sM"||put(SSNum,best.)),put(m,best.));
	call symputx (compress("sDSF"||put(SSNum,best.)),dsf);
	call symputx (compress("sNum"||put(SSNum,best.)),CRule);
	call symputx (compress("sCrt"||put(SSNum,best.)),CriteriaList);
run;
proc sql noprint;
	select count(dsf) into :sN trimmed
	from FinalSSCI;
quit;

* Create macro call for covariate var names 								;
*	&criteriaparen= "criteria1" "criteria2"... "criteriaN"					;
data temp_covarsa;run;
data temp_covarsa;
	length tempcovs1 tempcovs2 $500.;
	set temp_covarsa;

	tempcovs1="&covars";
	tempcovs2=left(trim('"'||strip(TRANWRD(trim(tempcovs1),trim(' '),'" "'))||'"'));
	call symput("covsparen",strip(tempcovs2));
run;
proc contents data=dsf.&dsn memtype=data 
	out=CovVars (keep=NAME LABEL VARNUM where=(NAME in (&covsparen))) 
	noprint;
run;
proc sort data=CovVars; by VARNUM;run;
data _null_;
	set CovVars;
	call symputx(compress("Cov"||put(_n_,best.)), NAME);
run;
proc sql noprint;
	select count(NAME) into :nCov trimmed
	from CovVars;
quit;
proc sort data=CovVars; by NAME;run;

*Run through the final DSFs from SS;
data DFCData;run;

%do i=1 %to &sN;

	*create dsf vars in nesarc;
	data &dsn._LR;
		set dsf.&dsn;
		if sum(of &&sCrt&i)>=&&sM&i then &&sDSF&i=1;
		else if .< sum(of &&sCrt&i) < &&sM&i then &&sDSF&i=0;
		*only disordered;
		if &dxvar=1 then output;
	run;

	data CovVarsX; set CovVars; run;
	*Screening Check to see if DFC should be assessed;
	%do j=1 %to &nCov;
		proc sort data=&dsn._LR; by &&Cov&j;run;
		proc surveyfreq data=&dsn._LR;
			weight &wtvar;
			by &&Cov&j;
			tables &&sDSF&i;
			ods output oneway=testss_Covars0 (keep=&&Cov&j &&sDSF&i Percent where=(&&sDSF&i=1)) ;
		run;
		data testss_Covars1;
			length NAME $32.;
			set testss_Covars0;
			NAME="&&Cov&j";
			if Percent >=97.0 then NoCheck=1;
		run;
		proc sql;
			create table testss_Covars2 as
			select NAME, count(&&Cov&j) as totcats, sum(NoCheck) as totNoCheck, 
				case 
					when calculated totNoCheck = calculated totcats then 1
					when calculated totNoCheck < calculated totcats then 0
					else .
				end as tmpNoCheck2
			from testss_Covars1
			group by NAME;
		quit;
		data CovVarsX;
			merge	CovVarsX (in=a)
						testss_Covars2 (in=b keep=NAME tmpNoCheck2);
			by NAME;
			if a;
			if b then do;
				NoCheck2=tmpNoCheck2;
			end;
			drop tmpNoCheck2;
		run;
		proc datasets lib = work;
			delete testss_:;
		quit;
		run;
	%end;

	*Determine if ANY Covars need to be checked;
	proc sql;
		select count(NoCheck2) into :totNoCheck trimmed
		from CovVarsX
		where NoCheck2=1;
	quit;

	%if "&totNoCheck"="&nCov" %then %do; *NO COVARS TO BE CHECKED - KEEP THE RULE;
		data temp_out3;run;
		data temp_out3;
			length dsf $60.;
			set 	temp_out3;
			dsf="&&sDSF&i";
			KeepDFC=1;
			CRule=&&sNum&i;
			p=&&sP&i;
			m=&&sM&i;
		run;
		data DFCData;
			set	DFCData
					temp_out3;
			if dsf = "" then delete;
		run;
	%end;

	%else %if "&totNoCheck" < "&nCov"%then %do; *THERE IS AT LEAST 1 COVAR TO CHECK DIFF;

		*Determine which covars to be assessed;
		*Macro vars for which covars to be used in LR;
		data CovVarsXx;
			length testvars $100.;
			retain testvars;
			set CovVarsX (where=(NoCheck2=0));
			by NAME;
			testvars=catx(" ", testvars, NAME);
			call symputx("testvars",testvars);
		run;
		ods html close;
		ods preferences;
		ods html newfile=none;
		
		*Proc SurveryLogistic;
		proc surveylogistic data=&dsn._LR;
			weight &wtvar;
			class &dxvar &covars;	
			model &&sDSF&i (event='1') = &dxvar &testvars;
			ods output ParameterEstimates=temp_out0;
		run;
		
		*determine if any covariates are significant;
		proc sql noprint;
			select count(ProbChiSq) into :nonsig trimmed
			from temp_out0
			where Variable ne "Intercept" and ProbChiSq >0.05;
		
			select count(ProbChiSq) into :totcov trimmed
			from temp_out0
			where Variable ne "Intercept";
		quit;
		%put &nonsig &totcov;
		
		data temp_out1;
			length Value dsf $60.;
			set temp_out0;
			Value=trim(left(compbl(put(Estimate,8.3)||"  /  "||left(trim(put(StdErr,8.3)))||"  /  "||left(trim(put(WaldChiSq,8.3)))||"  /  "||left(trim(put(ProbChiSq,8.3))))));
			dsf="&&sDSF&i";
			%if "&nonsig"="&totcov" %then %do;
				KeepDFC=1;
			%end;
			%else %if "&nonsig" ne "&totcov" %then %do;
				KeepDFC=0;
			%end;
		run;
		proc transpose data=temp_out1 out=temp_out2 (drop=_NAME_);
			by dsf KeepDFC;
			id Variable;
			var Value;
		run;
		data temp_out3;
			set temp_out2;
			CRule=&&sNum&i;
			p=&&sP&i;
			m=&&sM&i;
		run;
		data DFCData;
			set 	DFCData
					temp_out3 (drop=Intercept);
			if dsf = "" then delete;
		run;
	%end;
%end;

*DELETE DS;
proc datasets lib = work;
	delete temp:;
quit;
run;

*OUTPUT PERM;
data 	dsf.DFCData
		dsf.FinalDFC;
	set DFCData;
	output dsf.DFCData;
	if KeepDFC=1 then output dsf.FinalDFC;
run;

*FINAL DATA;

*turning log off;
options nonotes nomprint nosource nosource2 nosymbolgen nomlogic nomacrogen nosymbolgen nomlogic nomprint nomfile;
%**Write Message to Log***;
options notes;
%put;
%put NOTE: Finalize final list of optimal DSFs.;
%put;
options nonotes;

proc sort data=dsf.FinalCorrCI out=FinalCorrCI; by CRule;run;
proc sort data=dsf.FinalDFC out=FinalDFC0; by CRule p m;run;
proc sort data=dsf.FinalCombo out=FinalCombo; by CRule p;run;

data FinalDFC1;
	length corr95 $25.;
	merge	FinalDFC0 (in=a)
				FinalCorrCI (in=b keep=CRule PCorr:)
				FinalCombo (in=c keep=CRule CriteriaList CriteriaNumList C:);
	by CRule;
	if a;
	corr95 = compbl(trim(left(put(round(PCorr11Mean, 0.0001), 8.4)))||" ("||trim(left(put(round(PCorr11CI95L, 0.0001), 8.4)))||", "||
					trim(left(put(round(PCorr11CI95U, 0.0001), 8.4)))||")");
run;

proc sort data=FinalDFC1; by dsf;run;
proc sort data=dsf.FinalPrevCI out=FinalPrevCI; by dsf;run;
proc sort data=dsf.FinalSSCI out=FinalSSCI; by dsf;run;

data FinalDSF0;
	length sens95 spec95 prev95 prevDiff95 $25.;
	merge	FinalDFC1 (in=a)
				FinalPrevCI (in=b keep=dsf mPrev: stdPrev mPDiff: stdPDiff mPDiffAbs)
				FinalSSCI (in=c keep=dsf mSens: stdSens mSpec: stdSpec);
	by dsf;
	if a;

	SensSpec=mSens+mSpec;
	sens95 = compbl(trim(left(put(round(mSens, 0.0001), 8.4)))||" ("||trim(left(put(round(mSens95LB, 0.0001), 8.4)))||", "||trim(left(put(round(mSens95UB, 0.0001), 8.4)))||")");
	spec95 = compbl(trim(left(put(round(mSpec, 0.0001), 8.4)))||" ("||trim(left(put(round(mSpec95LB, 0.0001), 8.4)))||", "||trim(left(put(round(mSpec95UB, 0.0001), 8.4)))||")");
	prev95 = compbl(trim(left(put(round(mPrev, 0.0001), 8.4)))||" ("||trim(left(put(round(mPrev95L, 0.0001), 8.4)))||", "||trim(left(put(round(mPrev95U, 0.0001), 8.4)))||")");
	prevDiff95 = compbl(trim(left(put(round(mPDiff, 0.0001), 8.4)))||" ("||trim(left(put(round(mPDiff95L, 0.0001), 8.4)))||", "||trim(left(put(round(mPDiff95U, 0.0001), 8.4)))||")");
run;

proc sort data=FinalDSF0; by descending SensSpec descending PCorr11Mean;run;

data FinalDSF;
	set FinalDSF0;
	DSFn+1;
run;

data FinalVars;
	length finalvars $100.;
	retain finalvars;
	set CovVarsX;
	by NAME;
	finalvars=catx(",", finalvars, NAME);
	call symputx("finalvars",finalvars);
run;
proc sql;
	create table dsf.FinalDSF as
	select DSFn, dsf, p, m, CRule, corr95, prev95, sens95, spec95, keepDFC, &finalvars , CriteriaNumList, CriteriaList, 
		prevDiff95, PCorr11Mean, PCorr11CI95L, PCorr11CI95U, mPrev, stdPrev, mPrev95L, mPrev95U, mPDiff, stdPDiff, mPDiff95L, mPDiff95U, mPDiffAbs,
		mSens, stdSens, mSens95LB, mSens95UB, mSpec, stdSpec, mSpec95LB, mSpec95UB, mSpec1, SensSpec, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11
	from FinalDSF;
quit;

%**Write Message to Log***;
options notes;
%put;
%put NOTE: Macro Completed.;
%put;
options nonotes;

%end; *End DIFF Section;

%mend DSF;

