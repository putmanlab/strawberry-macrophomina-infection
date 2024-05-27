* ********************************** *
* STRAWBERRY Macrophomina infection  *
* Greenhouse						 *
* Pathogen quantification - plating	 *
* Roots - per sample (plant)		 *
* ********************************** *;

*** set variables for folder locations and base filename used for saving output
	** local to lab desktop;
*	%LET base_path=C:\Users\Alex Putman\GoogleDrive\UCR_VSPLab\07_SB_SCRI-pathology\SAS_Mac\infection;
*	%LET results_path=&base_path.\4_results\plating-samp-roots;
*	%LET results_path_img=&base_path.\4_results\plating-samp-roots\html_images;

	** to SAS ODA/SAS Studio;
	%LET base_path=/home/u63629950/sb-mac-inf;
	%LET results_path=&base_path./4_results/plating-samp-roots;
	%LET results_path_img=&base_path./4_results/plating-samp-roots/html_images;

	** both;
	%LET name_base=plating-roots_; 

*** load macros for controlling output;
	** local;
*	%include "&base_path.\3_analysis\output_files_macro.sas";
	
	** SAS Studio;
	%include "&base_path./3_analysis/output_files_macro.sas";

*options ls=120 nonumber formdlim=' ' pagesize=52 center;


* ********** *
* A. Import  *
* ********** *;

*** save log to file;
proc printto new log="&results_path/&name_base.A_sas-log.log"; run; 

** import data;
	* local;
*	proc import 
*			datafile="&base_path.\2_data\sb-mac-inf_data_sample.csv"
*			dbms=dlm replace out=plating_samp;
*		delimiter=",";
*		getnames=YES;
*	run;
			
	* SAS ODA/SAS Studio;
	proc import 
			datafile="&base_path./2_data/sb-mac-inf_data_sample.csv"
			dbms=dlm replace out=plating_samp;
		delimiter=",";
		getnames=YES;
	run;

** split data;
data samp_roots;
	set plating_samp;
	if tissue in ('crown') then delete;
run;

** sort dataset;
proc sort data=samp_roots;
	by trial treatment timepoint plant;
run;

** check dataset;
proc print data=samp_roots;
	title 'mac inf plating samp roots full review';
run;

*** direct log back to SAS Log window;
proc printto; run; 

** NOTE: code blocks used for initial testing and optimization were moved to appendices
	* whereas code used in the final analysis workflow (e.g., that random effect term was selected) were not moved but commented out because subsequent steps provide same results


* ******************************** *
* B. Analyze - Identify Base Model *
* ******************************** *;

*** Step 1 ***;
	** OBJ: simple model, test distributions;

	** Step 1-1 **;
		* OBJ: negative binomial;
		* see Appendix B for code
		* RESULTS: all pearson/DF < 1
			* 1: treatment, treatment*timepoint
			* 2: none signif
			* 3: weak treatment*timepoint
			* 4: treatment*timepoint

	** Step 1-2 **;
		* OBJ: poisson;
		* see Appendix B for code
		* RESULTS: all pearson/DF very high (>50)
		
	** CONCLUSION: does not make sense, cannot detect difference with control which are all 0s
		* adjust values to remove 0s


* ******************************************* *
* C. Analyze - Adjusted - Identify Base Model *
* ******************************************* *;

*** adjust values;
	** previous rounds of analyses failed with raw data, add 1 to all values to remove 0s;
data samp_roots_adj;
	set samp_roots;
	colonies_sum_per_g_adj = colonies_sum_per_g + 1;
	drop colonies_sum_per_g;
run;

*** Step 1 ***;
	** OBJ: simple model, test distributions;

	** Step 1-1 **;
		* OBJ: negative binomial;
			* NOTE: commented out because later section provides the same results;
*		%output_html_files_on(name_step=C_step-1-1, title_step=C step 1-1); * redirect html, log output to files;
*		proc glimmix data=samp_roots_adj plot=residualpanel method=laplace;
*			class trial treatment timepoint;
*			by trial;
*			model colonies_sum_per_g_adj = treatment|timepoint / dist=negbinomial link=log htype=3;
*			title 'mac inf plating samp roots C adj id base model - step 1 simple, distributions - step 1-1 negative binomial';
*		run;
*		%output_html_files_off(); * turn off redirecting html, log output to files;

		* RESULTS: all pearson/DF < 1
			* treatment*timepoint signif for all four trials

	** Step 1-2 **;
		* OBJ: poisson;
		* see Appendix C for code
		* RESULTS: all pearson/DF very high (>50)
		
	** CONCLUSION: use adjusted values, negative binomial distribution


* ************************************ *
* D. Analyze - Adjusted - Main effects *
* ************************************ *;

*** Step 1 ***;
	** OBJ: analyze fixed effects with method=laplace removed and ddfm=kr2 (from B step 3-1);

	** Step 1-1 **;
		* OBJ: analyze;
			* NOTE: commented out because later section provides the same results;
*		%output_html_files_on(name_step=D_step-1-1, title_step=D step 1-1); * redirect html, log output to files;
*		proc glimmix data=samp_roots_adj plot=residualpanel;
*			class trial treatment timepoint;
*			by trial;
*			model colonies_sum_per_g_adj = treatment|timepoint / dist=negbinomial link=log htype=3 ddfm=kr2;
*			title 'mac inf plating samp roots D adj main effect - step 1 -lapace, +kr2 - step 1-1 analyze';
*		run;
*		%output_html_files_off(); * turn off redirecting html, log output to files;

	** RESULTS: treatment*timepoint signif for all four trials
	** CONCLUSION: investigate interaction


* *********************************** *
* E. Analyze - Adjusted - Interaction *
* *********************************** *;

*** Step 1 ***;
	** OBJ: investigate interaction;

	** Step 1-1 **;
		* OBJ: investigate;
			* NOTE: commented out because later section provides the same results;
*		%output_html_files_on(name_step=E_step-1-1, title_step=E step 1-1); * redirect html, log output to files;
*		proc glimmix data=samp_roots_adj;
*			class trial treatment timepoint;
*			by trial;
*			model colonies_sum_per_g_adj = treatment|timepoint / dist=negbinomial link=log htype=3 ddfm=kr2;
*			
*			slice treatment*timepoint / sliceBy=treatment;
*			
*			title 'mac inf plating samp roots E adj interaction - step 1 investigate - step 1-1 investigate';
*			ods exclude MeanPlot SliceDiffs DiffPlot;
*		run;
*		%output_html_files_off(); * turn off redirecting html, log output to files;

	** RESULTS: effect of timepoint significant within crown and root inoculation treatments for all trials
	** CONCLUSION: separate means


* ************************************* *
* F. Analyze - Adjusted - Examine means *
* ************************************* *;

*** Step 1 ***;
	** OBJ: examine means;

	** Step 1-1 **;
		* OBJ: examine;
		%output_html_files_on(name_step=F_step-1-1, title_step=F step 1-1); * redirect html, log output to files;
		proc glimmix data=samp_roots_adj;
			class trial treatment timepoint;
			by trial;
			model colonies_sum_per_g_adj = treatment|timepoint / dist=negbinomial link=log htype=3 ddfm=kr2;
			
			slice treatment*timepoint / sliceBy=treatment means ilink linestable adjust=tukey;
			
			title 'mac inf plating samp roots F adj means - step 1 examine - step 1-1 examine';
			ods exclude MeanPlot SliceDiffs DiffPlot;
		run;
		%output_html_files_off(); * turn off redirecting html, log output to files;

	** CONCLUSION: analysis complete



* --------------------APPENDICES-------------------- *;



* ***************************************** *
* Appendix B. Analyze - Identify Base Model *
* ***************************************** *;

*** Step 1 ***;
	** OBJ: simple model, test distributions;

	** Step 1-1 **;
		* OBJ: negative binomial;
*		%output_html_files_on(name_step=B_step-1-1, title_step=B step 1-1); * redirect html, log output to files;
*		proc glimmix data=samp_roots plot=residualpanel method=laplace;
*			class trial treatment timepoint;
*			by trial;
*			model colonies_sum_per_g = treatment|timepoint / dist=negbinomial link=log htype=3;
*			title 'mac inf plating samp roots B id base model - step 1 simple, distributions - step 1-1 negative binomial';
*		run;
*		%output_html_files_off(); * turn off redirecting html, log output to files;
**
*	** Step 1-2 **;
*		* OBJ: poisson;
*		%output_html_files_on(name_step=B_step-1-2, title_step=B step 1-2); * redirect html, log output to files;
*		proc glimmix data=samp_roots plot=residualpanel method=laplace;
*			class trial treatment timepoint;
*			by trial;
*			model colonies_sum_per_g = treatment|timepoint / dist=poisson link=log htype=3;
*			title 'mac inf plating samp roots B id base model - step 1 simple, distributions - step 1-2 poisson';
*		run;
*		%output_html_files_off(); * turn off redirecting html, log output to files;


* **************************************************** *
* Appendix C. Analyze - Adjusted - Identify Base Model *
* **************************************************** *;

*** Step 1 ***;
	** OBJ: simple model, test distributions;

	** Step 1-2 **;
		* OBJ: poisson;
*		%output_html_files_on(name_step=C_step-1-2, title_step=C step 1-2); * redirect html, log output to files;
*		proc glimmix data=samp_roots_adj plot=residualpanel method=laplace;
*			class trial treatment timepoint;
*			by trial;
*			model colonies_sum_per_g_adj = treatment|timepoint / dist=poisson link=log htype=3;
*			title 'mac inf plating samp roots C adj id base model - step 1 simple, distributions - step 1-2 poisson';
*		run;
*		%output_html_files_off(); * turn off redirecting html, log output to files;

	
