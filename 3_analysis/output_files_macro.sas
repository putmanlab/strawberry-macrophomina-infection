* ********************************** *
* Putman lab UC Riverside			 *
* SAS macros						 *
* saving output to files			 *
* ********************************** *;


*** save html output and log as separate files for each step;
	** turn on (used before proc call);
	%macro output_html_files_on(name_step, title_step);
		* redirect html output to file, set html file names and paths;
		ods html body="&name_base.&name_step._sas-output.html" (title="&title_step.") path="&results_path" gpath="&results_path_img" (url="html_images/"); 
		* turn html graphics on, have html image suffixes reset after each name_step;
		ods graphics on / reset imagename="&name_step.";
		* redirect log output to file, set filename;
		proc printto new log="&results_path/&name_base.&name_step._sas-log.log"; run;
	%mend;

	** turn off saving output for that step;
	%macro output_html_files_off();
		* turn graphics off;
		ods graphics off;
		* close html;
		ods html close;
		* redirect log output to regular window;
		proc printto; run;
	%mend;	

*** control regular proc output;
	** turn off (used before proc call);
	%macro output_off();
		ods graphics off; 
		ods exclude all; 
		ods noresults;
	%mend;
	
	** turn on (used after proc call);
	%macro output_on();
		ods graphics on; 
		ods exclude none; 
		ods results;
	%mend;

*** export combined output tables to csv;
	%macro output_table_csv(table_obj, name_step);
		proc export data=&table_obj.
			outfile="&results_path./&name_base.&name_step._sas-output_table.csv"
			dbms=csv 
			replace;
		run;
	%mend;