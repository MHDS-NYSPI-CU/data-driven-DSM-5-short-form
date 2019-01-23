
%macro DSF_TF_v2	(	file=,	/* REQUIRED - File location for all datasets created using Macro DSF_v2 */
									dsn=, /* REQUIRED - name of original pariticipant level dataset */
									dxvar=,	/* REQUIRED - dsm diagnosis variable */
									sudlab=,/*OPTIONAL - label for histogram (i.e., sudlab=MAUD)*/
									wtvar=,	/* OPTIONAL - name of weight variable in original participant level dataset - default=(missing) - no weight used unless variable is identified */
									p1cut=,	/*OPTIONAL - Prevalence Plot 1 cut point for DSFs (DEFAULT = EMPLY --? ALL Shown - can put a number 200 so only the first 200 DSFs will be displayed */
									adjlabel=, /* OPTIONAL - Final DSF figure for prevalence - this adjusts the prevalence label to the right (+ #) or left (- #) depending on how you need it */
									size=medium, /*OPTIONAL -  if there are many DSFs than choose size=large, if a medium amount then choose size=medium, if low number of DSFs choose size=small (DEFAULT=medium)*/
									crtlabel= /* OPTIONAL - if you want to add criterion labels to criteria prevalence histogram then crtlabel=[name of format you decide on] and must create labels in proc format provided */
);

/*FORMATS*/
proc format;
value criteriaf
	0="Excluded"
	1="Included"
;
run;
options ls=80 ps=55 nodate nofmterr mprint mlogic symbolgen MergeNoBy=error MSGLevel = I compress=yes;
libname dsf "&file";

*******************************;
*% of current users that have disorder (95% CIs);
*******************************;
proc means data=dsf.BootIDData1 noprint;
	%if %length(&wtvar)>0 %then %do;
		weight &wtvar;
	%end;
	by BootSamp;
	vars &dxvar;
	output out=Prev0 mean=Prev;
run;
proc means data=Prev0;
	vars Prev;
	output out=Prev1 mean=Prev std=std;
run;
data Prev;
	set Prev1;

	*95% CI;
	PrevCI95L = Prev - (1.96*std);
	PrevCI95U = Prev + (1.96*std);
run;

*******************************;
* Histogram for Criteria Prevalence in Original Sample;
*******************************;

*obtain macro call for list of criteria;
data CriteriaVars0;
	set dsf.CriteriaVars;
	length=length(NAME);
run;
proc sort data=CriteriaVars0; by length NAME;run;
data CriteriaVars;
	length CriteriaVars $100.;
	retain CriteriaVars;
	set CriteriaVars0;
	call symputx(compress("C"||put(_n_,best.)), NAME);
	by length NAME;
	CriteriaVars=catx(" ", CriteriaVars, NAME);
	call symputx("CriteriaVars",CriteriaVars);
run;
proc sql noprint;
	select count(NAME) into :totCrt trimmed
	from CriteriaVars0;
quit;

*Prevelance of criteria in original sample;
proc means data=dsf.&dsn n mean;
	%if %length(&wtvar)>0 %then %do;
		weight &wtvar;
	%end;
	var &CriteriaVars;;
	output out=prevalence0;
run;

*Create data for output;
data prevalence1 (keep = Criterion Prevalence prevc);
	length prevc $15.;
	set prevalence0 (where = (_STAT_ = "MEAN"));
	%macro prev;
		%do i=1 %to &totCrt;
			Criterion=&i;
			Prevalence=round((&&C&i)*100,0.1);
			prevc=left(trim(put(Prevalence,8.1)||"%"));
			output;
		%end;
	%mend;
	%prev;
run;

ods html close;
ods preferences;
ods html newfile=none;

proc print data=Prev;
title1 "prevalence of diagnosis among current users with 95% CIs";
run;

proc freq data=dsf.&dsn;
	table &dxvar;
title1 "prevalence of diagnosis among current users (raw, unweighted data - to determine raw totals)";
run;
title;

ods graphics / height=5.5in width=8in;
proc sgplot data= Prevalence1;
	vbar Criterion / response = Prevalence datalabel = prevc;
	xaxis label = "DSM-5 Criteria" discreteorder=data;
	yaxis label = "Prevalence (%) of Criteria Among Current Users";
	%if %length(&crtlabel)>0 %then %do;
		format Criterion &crtlabel..;
	%end;
run;
ods graphics off;

*******************************;
* Boxplots of Weighted Pearson Correlations for all Subscales (figure 3.5);
*******************************;

ods graphics on;
ods graphics / height=5in width=7in;
proc sgplot data=dsf.CorrCIData;
	vbox PCorr11Mean / category=p ;
	xaxis label = "Number of Criteria in Subscale (p)" grid;
	yaxis label = "Pearson Correlations with Total" grid;
run;
ods graphics off;
title;

*******************************;
* PREV PLOT 1&2 - Scatterplot of Prevalences (figure 4.6) & Panelby p for the difference (AUD - prev) and 95% CIs for the top chosen within 5% prev (figure 4.7);
*******************************;
proc template;
	define style prev1;
		Parent=styles.htmlblue;
		style Graph from Graph / attrpriority="None";
		style GraphBar from GraphComponent / displayopts="fill outline fillpattern";
		style GraphData1 from GraphData1 / color=blue markersymbol="circle";
		style GraphData2 from GraphData2 / color=red markersymbol="star";
		style GraphData3 from GraphData3 / color=purple markersymbol="square";
		style GraphData4 from GraphData4 / color=orange markersymbol="starfilled";
		style GraphData5 from GraphData5 / color=grey markersymbol="diamond";
		style GraphData6 from GraphData6 / color=gold markersymbol="plus";
		style GraphData7 from GraphData7 / color=yellow markersymbol="circlefilled";
		style GraphData8 from GraphData8 / color=magenta markersymbol="hash";
		style GraphData9 from GraphData9 / color=green markersymbol="triangle";
		style GraphData10 from GraphData10 / color=cyan markersymbol="trianglefilled";
	end;
run;

*********;
*PREV PLOT 1;
*********;
*DSM-5 Reference Values & 0.05 Reference Lines!;
proc sql noprint;
	select mDSMPrev, mDSMPrevCI95LB, mDSMPrevCI95UB, sum(mDSMPrev+0.05), sum(mDSMPrev-0.05), sum(mDSMPrev &adjlabel), round(mDSMPrev*100,0.1)
		into :DSMp trimmed, :DSMlb trimmed, :DSMub trimmed, :ub5 trimmed, :lb5 trimmed, :prevref trimmed, :prevperc trimmed
	from dsf.DSMPrev;
quit;

*Creat Temp Rank var;
proc sort data=dsf.PrevCIData out=PrevCIData0; by mPDiffAbs;run;
data PrevCIData;
	set PrevCIData0;
	by mPDiffAbs;
	Rank+1;
run;
proc sort data=PrevCIData; by p;run;
ods html style=prev1;
ods graphics on / height=5in width=8in;
proc sgplot data=PrevCIData;
%if %length(&p1cut)>0 %then %do;
	where Rank <=&p1cut;
%end;
	scatter y=mPrev 	x=Rank / 	group=p name="main" 
													errorbarattrs=(thickness=1pt pattern=dash) yerrorlower=mPrev95L yerrorupper=mPrev95U;
	refline &DSMp &DSMlb &DSMub / lineattrs=(color=gray thickness=1pt pattern=dash) name='ref1' ;
	refline &ub5 &lb5/ lineattrs=(color=red thickness=1pt pattern=dash) name='ref2' ;

	keylegend "main" /title="Short-Form Size (p)";
	%if %length(&sudlab)=0 %then %do;
		keylegend "ref1" / title="DSM-5 Disorder Prevalence (95% CIs)" location=inside position=topleft across=1;
		keylegend "ref2" / title="Within 5% of DSM-5 Disorder Prevalence" location=inside position=bottomleft across=1;
	%end;
	%if %length(&sudlab)>0 %then %do;
		keylegend "ref1" / title="DSM-5 &sudlab Prevalence (95% CIs)" location=inside position=topleft across=1;
		keylegend "ref2" / title="Within 5% of DSM-5 &sudlab Prevalence" location=inside position=bottomleft across=1;
	%end;
	xaxis label = "Rank of Diagnostic Short-Forms" grid;
	yaxis label = "DSF Prevalence" grid;
run;
ods graphics off;
title;

*********;
*PREV PLOT 2;
*********;
proc sort data=dsf.FinalPrevCI out=FinalPrevCI; by mPDiffAbs;run;
ods html style=styles.htmlblue;
ods graphics on / height=5in width=8in;
proc sgpanel data=FinalPrevCI;
	panelby p / columns=10 colheaderpos=top uniscale=row;
	scatter y=mPDiff 	x=CRule /	group=m name="main" 
													errorbarattrs=(thickness=1pt pattern=dash) yerrorlower=mPDiff95L yerrorupper=mPDiff95U;
	refline -0.05 0.05/ lineattrs=(color=red thickness=1pt pattern=dash) name='ref1' ;
	keylegend "main" /title="Threshold (m)";
	%if %length(&sudlab)=0 %then %do;
		keylegend "ref1" / title="Within 5% of DSM-5 Disorder Prevalence" position=bottom ;
	%end;
	%if %length(&sudlab)>0 %then %do;
		keylegend "ref1" / title="Within 5% of DSM-5 &sudlab Prevalence" position=bottom ;
	%end;
	colaxis grid display=none;
	rowaxis label = "Difference in Prevalence" grid;
run;
ods graphics off;
title;

*******************************;
* SENS/SPEC - (figure 4.8);
*******************************;
proc sort data=dsf.SSCIData out=SSCIData; by p;run;

ods html style=prev1;
ods graphics on / height=5in width=6in;
proc sgplot data=SSCIData;
	scatter y=mSens 	x=mSpec1 / 	group=p name="main" 
																errorbarattrs=(thickness=1pt pattern=dash) yerrorlower=mSens95LB yerrorupper=mSens95UB;
	keylegend "main" /title="Number of Criteria (p)";
	xaxis label = "1- Specificity (FPR) = Pr(DSF + | DSM5 - )" values=(0 to 0.16 by 0.02) grid;
	yaxis label = "Sensitivity (TPR) = Pr(DSF + | DSM5 +)" grid values=(0.65 to 1 by 0.05) grid;
run;
ods graphics off;
title;

*******************************;
* HEATMAP OF FINAL DSF;
*******************************;
proc sort data=dsf.FinalDSF out=FinalDSF; by descending SensSpec;run;

data FinalDSF_Long1;
	length Criteria $20.;
	set FinalDSF;
	if CRule ne . then do;
		%macro long;
			%do i=1 %to &totCrt;
				if C&i=0 then do;
					Criteria="C&i";
					sort=&i;
					Criteria2=0;
					output;
				end;
				else if C&i=1 then do;
					Criteria="C&i";
					sort=&i;
					Criteria2=1;
					output;
				end;
			%end;
		%mend;
		%long;
	end;
run;

*Finalize Long data for heatmap;
data FinalDSF_Long;
	length DSFnx temp1 temp2 sortx $50.;
	set FinalDSF_Long1;
	temp1=trim(left(compbl("p="||compress(left(put(p,8.))||","))));
	temp2=trim(left(compbl("m="||compress(left(put(m,8.))))));
	DSFnx=compbl(trim(left("DSF "||compbl(trim(left(put(CRule,8.)))||": ")||temp1||temp2)));

	sortx=trim(put(sort,labelf.));

	%if %length(&crtlabel)>0 %then %do;
		format Criteria2 criteriaf. sort sortcfx.;
	%end;
	%if %length(&crtlabel)=0 %then %do;
		format Criteria2 criteriaf. ;
	%end;
	%if %length(&sudlab)=0 %then %do;
		label sortx = "DSM-5 Criteria";
	%end;
	%if %length(&sudlab)>0 %then %do;
		label sortx = "DSM-5 &sudlab Criteria";
	%end;
run;

*CHERI - I had to create the sort option so I could order the criteria least --> most prev in heat map;
*HeatMap Template;
proc template;
	define statgraph heat;
		begingraph;
			layout overlay / 
					xaxisopts=(type=discrete labelattrs=(size=10pt) discreteopts=(tickvaluefitpolicy=rotate) tickvalueattrs=(size=10pt))
					%if "&size"="large" %then %do;
						yaxisopts=(label="Optimal DSFs (ordered by Sensitivity/Specificity)" labelattrs=(size=10pt) tickvalueattrs=(size=7.5pt) reverse=true);
					%end;
					%else %do;
						yaxisopts=(label="Optimal DSFs (ordered by Sensitivity/Specificity)" labelattrs=(size=10pt) tickvalueattrs=(size=9pt) reverse=true);
					%end;
				heatmapparm y=DSFnx x=sortx colorgroup=criteria2 / 
					name='heat' 
					display=all /*Displays grid around filled blocks*/
					fillattrs=(transparency=.2) /*transparency of black (on)*/
					;
				discretelegend 'heat' / location=outside valign=bottom title="Criteria Included/Excluded" valueattrs=(size=10pt);
			endlayout;
		endgraph;
	end;

	define style colors1;
		parent=styles.htmlblue;
		class graphdata1 / color=black;
		class graphdata2 / color=lightgray;
	end;
	define style colors2;
		parent=styles.htmlblue;
		class graphdata1 / color=lightgray;
		class graphdata2 / color=black;
	end;
run;

ods graphics on;
ods html style=colors1;
%if "&size"="large" %then %do;
	ods graphics / reset antialiasmax=7000 width=8in height=8in;
%end;
%if "&size"="medium" %then %do;
	ods graphics / reset antialiasmax=7000 width=8in height=6.5in;
%end;
%if "&size"="small" %then %do;
	ods graphics / reset antialiasmax=7000 width=8in height=5in;
%end;
proc sgrender data=FinalDSF_Long template=heat;
run;

*******************************;
* TABLE OF FINAL FIT STATS;
*******************************;
proc print data=FinalDSF;
	var CRule p m corr95 prev95 sens95 spec95;
run;

*******************************;
* FIGURES OF FINAL FIT STATS;
*******************************;
data finallab;
	length label $25.;
	set FinalDSF;
	label=left(trim(compbl("DSF: "||compress(left(trim(put(CRule,best.)))||",")||" p="||compress(left(trim(put(p,best.)))||",")||" m="||left(trim(put(m,best.))))));
run;
proc sort data=finallab; by SensSpec;run;
%if "&size"="large" %then %do;
	ods graphics on / height=6.5in width=4in;
%end;
%if "&size"="medium" %then %do;
	ods graphics on / height=4.5in width=4in;
%end;
%if "&size"="small" %then %do;
	ods graphics on / height=3.5in width=4in;
%end;
proc sgplot data=finallab ;
	scatter y=label 	x=PCorr11Mean / 	name="main" 
														errorbarattrs=(thickness=1pt pattern=dash) xerrorlower=PCorr11CI95L xerrorupper=PCorr11CI95U;

	xaxis label = "Correlation with Total" grid reverse ;
	%if "&size"="large" %then %do;
		%if %length(&sudlab)=0 %then %do;
			yaxis label = "Optimal DSFs" grid discreteorder=data valueattrs=(size=7.5pt);
		%end;
		%if %length(&sudlab)>0 %then %do;
			yaxis label = "Optimal DSFs for DSM-5 &sudlab" grid discreteorder=data valueattrs=(size=7.5pt);
		%end;
	%end;
	%else %do;
		%if %length(&sudlab)=0 %then %do;
			yaxis label = "Optimal DSFs" grid discreteorder=data ;
		%end;
		%if %length(&sudlab)>0 %then %do;
			yaxis label = "Optimal DSFs for DSM-5 &sudlab" grid discreteorder=data;
		%end;
	%end;
run;
ods graphics off;
%if "&size"="large" %then %do;
	ods graphics on / height=6.5in width=2.5in;
%end;
%if "&size"="medium" %then %do;
	ods graphics on / height=4.5in width=2.5in;
%end;
%if "&size"="small" %then %do;
	ods graphics on / height=3.5in width=2.5in;
%end;
proc sgplot data=finallab ;
	scatter y=label 	x=mSens / 	name="main" 
														errorbarattrs=(thickness=1pt pattern=dash) xerrorlower=mSens95LB xerrorupper=mSens95UB;

	xaxis label = "Sensitivity" grid reverse ;
	yaxis display=none grid;
run;
ods graphics off;

%if "&size"="large" %then %do;
	ods graphics on / height=6.5in width=2.5in;
%end;
%if "&size"="medium" %then %do;
	ods graphics on / height=4.5in width=2.5in;
%end;
%if "&size"="small" %then %do;
	ods graphics on / height=3.5in width=2.5in;
%end;
proc sgplot data=finallab ;
	scatter y=label 	x=mSpec/ 	name="main" 
														errorbarattrs=(thickness=1pt pattern=dash) xerrorlower=mSpec95LB xerrorupper=mSpec95UB;

	xaxis label = "Specificity" grid reverse ;
	yaxis display=none grid;
run;
ods graphics off;

%if "&size"="large" %then %do;
	ods graphics on / height=6.5in width=2.5in;
%end;
%if "&size"="medium" %then %do;
	ods graphics on / height=4.5in width=2.5in;
%end;
%if "&size"="small" %then %do;
	ods graphics on / height=3.5in width=2.5in;
%end;
proc sgplot data=finallab ;
	scatter y=label 	x=mPrev/ 	name="main" 
														errorbarattrs=(thickness=1pt pattern=dash) xerrorlower=mPrev95L xerrorupper=mPrev95U;
	refline &DSMp / axis=x lineattrs=(color=black);
	%if %length(&sudlab)=0 %then %do;
		refline &prevref / axis=x label=("Prev (&prevperc.%)") labelloc=inside transparency=1;
	%end;
	%if %length(&sudlab)>0 %then %do;
		refline &prevref / axis=x label=("&sudlab Prev (&prevperc.%)") labelloc=inside transparency=1;
	%end;
	xaxis label = "Prevalence" grid ;
	yaxis display=none grid;
run;
ods graphics off;
%mend DSF_TF_v2;
