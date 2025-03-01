######################################################
# STRAWBERRY Macrophomina infection    			 	 #
# Greenhouse			 						 	 #
# Pathogen quantification - plating					 #
######################################################

## built on Docker putmanlab/exploratory-analysis:420.0

library(conflicted)

library(dplyr)
library(forcats)
library(ggplot2)
library(lubridate)
library(readr)
library(stringr)
library(tidyr)

library(patchwork) # for wrap_plots

conflict_prefer("date", "lubridate")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("spread", "tidyr")

setwd("/home/strawberry-macrophomina-infection")


################
#### Import ####
################

### upload
	## colonies
	colonies1 = read_csv(file="./2_data/October 2021 Layered Inoculum Study - Colonies.csv", col_names=T, na=c("NA"))
	colonies2 = read_csv(file="./2_data/November 2021 Layered Inoculum Study - Colonies.csv", col_names=T, na=c("NA"))
	colonies3 = read_csv(file="./2_data/February 2022 Layered Inoculum Study - Colonies.csv", col_names=T, na=c("NA"))
	colonies4 = read_csv(file="./2_data/March 2022 Layered Inoculum Study  - Colonies.csv", col_names=T, na=c("NA"))

	## tissue
	tissue1 = read_csv(file="./2_data/October 2021 Layered Inoculum Study - Tissue.csv", col_names=T, na=c("NA"))
	tissue2 = read_csv(file="./2_data/November 2021 Layered Inoculum Study - Tissue.csv", col_names=T, na=c("NA"))
	tissue3 = read_csv(file="./2_data/February 2022 Layered Inoculum Study - Tissue.csv", col_names=T, na=c("NA"))
	tissue4 = read_csv(file="./2_data/March 2022 Layered Inoculum Study  - Tissue.csv", col_names=T, na=c("NA"))

### combine trials and data
	## join trials
	in.1 = colonies1 %>% left_join(tissue1, by=c(c("Timepoint","Treatment","Tissue","Replicate","Plate")) )
	in.2 = colonies2 %>% left_join(tissue2, by=c(c("Timepoint","Treatment","Tissue","Replicate","Plate")) )
	in.3 = colonies3 %>% left_join(tissue3, by=c(c("Timepoint","Treatment","Tissue","Replicate","Plate")) )
	in.4 = colonies4 %>% left_join(tissue4, by=c(c("Timepoint","Treatment","Tissue","Replicate","Plate")) )

	## Add column to indicate trial
	in.1 = in.1 %>% mutate(trial="1")
	in.2 = in.2 %>% mutate(trial="2")
	in.3 = in.3 %>% mutate(trial="3")
	in.4 = in.4 %>% mutate(trial="4")
	
	## bind
	in.data = bind_rows(in.1, in.2, in.3, in.4)
	
### organize, clean up df
	## change column names
	in.data = in.data %>% rename(timepoint=Timepoint, treatment_old=Treatment, tissue_old=Tissue, plant=Replicate, plate=Plate, g_tissue_total=`Amount of Tissue (g)`, n_colonies=`# Mp Colonies`)

	## reduce treatment names, lower case tissue
		# make new column
		in.data = in.data %>% mutate(
			treatment=case_when(
				(treatment_old == "Control") ~ "control",
				(treatment_old == "Crown Inoculation") ~ "crown",
				(treatment_old == "Root Inoculation") ~ "root"),
			tissue=case_when(
				(tissue_old == "Crown") ~ "crown",
				(tissue_old == "Roots") ~ "roots") )
		
		# remove old column
		in.data = in.data %>% select(-treatment_old, -tissue_old)

	## change column types
	in.data = in.data %>% mutate(timepoint=as.integer(timepoint), n_colonies=as.integer(n_colonies) )
	
	## replace 0s; some samples showed a weight of 0 on scale, but minimum measurable weight is 0.01 g
	in.data = in.data %>% mutate(g_tissue_total=replace(g_tissue_total, g_tissue_total == 0, 0.009) )

### organize - per plate
	## replace g_tissue with per plate value; value entered is for whole sample, which was distributed across 7 plates
		# create new column
		data.inf.plate = in.data %>% mutate(g_tissue_plate=round(g_tissue_total / 7, digits=3) )
		
		# remove old column
		data.inf.plate = data.inf.plate %>% select(-g_tissue_total)
	
	## calculate colonies per g
	data.inf.plate = data.inf.plate %>% mutate(colonies_per_g=round(n_colonies / g_tissue_plate, digits=0) )

	## reorder columns
	data.inf.plate = data.inf.plate %>% select(trial, treatment, timepoint, plant, tissue, plate, g_tissue_plate, n_colonies, colonies_per_g)

### organize - per sample
	## checks
		# for missing plates
		in.data %>% group_by(trial, treatment, timepoint, plant, tissue) %>% summarize(ct=n(), ct_plate=n_distinct(plate) ) %>% filter(ct != 7 | is.na(ct) | ct_plate != 7 | is.na(ct_plate) )

		# one g_tissue_total per sample
		in.data %>% group_by(trial, treatment, timepoint, plant, tissue) %>% summarize(n_g_tissue=n_distinct(g_tissue_total) ) %>% filter(n_g_tissue > 1)
	
	## sum colonies for each sample
	data.inf.samp = in.data %>% 
		group_by(trial, treatment, timepoint, plant, tissue, g_tissue_total) %>%
		summarize(n_colonies_sum=sum(n_colonies) ) %>%
		ungroup()
		
	## calculate colonies per g
	data.inf.samp = data.inf.samp %>% mutate(colonies_sum_per_g=round(n_colonies_sum / g_tissue_total, digits=0) )
		
	## sort
	data.inf.samp = data.inf.samp %>% arrange(trial, treatment, timepoint, plant, tissue)
	
### export
	write_csv(data.inf.plate, file="./2_data_curated/sb-mac-inf_data_plate.csv", col_names=T, na=".", append=F)	
	write_csv(data.inf.samp, file="./2_data_curated/sb-mac-inf_data_sample.csv", col_names=T, na=".", append=F)	


###################
#### Summarize ####
###################
# summary statistics for Results

### 1. distribution of sample weights 
	## get total number of samples
		summ.1.tot = data.inf.samp %>% 
			distinct(trial, treatment, timepoint, plant, tissue) %>%
			group_by(trial, tissue) %>%
			summarize(n_total=n() )
		
		summ.1.tot.b = data.inf.samp %>% 
			distinct(trial, treatment, timepoint, plant, tissue) %>%
			group_by(tissue) %>%
			summarize(n_total=n() )
	
	## determine distribution
		# by each weight; collapse to plant samples, count number of samples for each weight, join totals, calculate percentage, spread for viewing
		summ.1 = data.inf.samp %>% 
			group_by(trial, tissue, g_tissue_total) %>%
			summarize(n_colonies=n()) %>%
			ungroup() 

		# by weight group
			# add grouping column
			summ.1.grp = summ.1 %>% mutate(weight_group=case_when(
				(g_tissue_total == 0.15) ~ "0.15",
				(g_tissue_total >= 0.10 & g_tissue_total < 0.15) ~ "0.10 to 0.14",
				(g_tissue_total >= 0.05 & g_tissue_total < 0.10) ~ "0.05 to 0.09",
				(g_tissue_total >= 0.00 & g_tissue_total < 0.05) ~ "0.00 to 0.04") )
						
	## export/show summaries
		# spread counts - supplementary table
		summ.1.ct = summ.1 %>%
			unite("tissue_trial", c("tissue","trial"), sep="_", remove=FALSE) %>%
			select(-tissue, -trial) %>%
			spread(key=tissue_trial, value=n_colonies) %>%
			mutate(across(c(crown_1:roots_4), ~ replace_na(.x, 0) ) ) 
			
			# export
			write_csv(summ.1.ct, file="./4_output/01_summary_tissue-weight-counts.csv", col_names=T, na="NA", append=F)	
		
		# weight groups
			# by trial - for results text (roots)
			summ.1.grp %>% 
				group_by(trial, tissue, weight_group) %>%
				summarize(n_colonies=sum(n_colonies) ) %>%
				ungroup() %>%
				left_join(summ.1.tot, by=c(c("trial","tissue")) ) %>%
				mutate(perc_colonies=round(n_colonies / n_total, digits=3) * 100) %>%
				unite("tissue_trial", c("tissue","trial"), sep="_", remove=FALSE) %>%
				select(-trial, -tissue, -n_colonies, -n_total) %>%
				spread(key=tissue_trial, value=perc_colonies) %>%
				mutate(across(c(crown_1:roots_4), ~ replace_na(.x, 0) ) ) 		
					
			# overall - for results text (crown)
			summ.1.grp %>% 
				group_by(tissue, weight_group) %>%
				summarize(n_colonies=sum(n_colonies) ) %>%
				ungroup() %>%
				left_join(summ.1.tot.b, by=c(c("tissue")) ) %>%
				mutate(perc_colonies=round(n_colonies / n_total, digits=3) * 100) %>%
				select(-n_colonies, -n_total) %>%
				spread(key=tissue, value=perc_colonies)
				
		# convert to percent and spread - for results text
		summ.1 %>%
			left_join(summ.1.tot, by=c(c("tissue","trial")) ) %>%
			mutate(percent=round(n_colonies / n_total, digits=3) * 100 ) %>%
			unite("tissue_trial", c("tissue","trial"), sep="_", remove=FALSE) %>%
			select(-tissue, -trial, -n_colonies, -n_total) %>%
			spread(key=tissue_trial, value=percent) %>%
			mutate(across(c(crown_1:roots_4), ~ replace_na(.x, 0) ) ) 		

### 2. by sample - median per timepoint - for results text
	data.inf.samp %>%
		group_by(trial, treatment, timepoint, tissue) %>%
		summarize(col_per_g_median=median(colonies_sum_per_g) ) %>%
		ungroup() %>%
		mutate(tissue_timepoint=paste(tissue, sprintf("%02d", timepoint), sep="_") ) %>%
		select(-timepoint, -tissue) %>%
		spread(key=tissue_timepoint, value=col_per_g_median) %>%
		arrange(treatment, trial)

### 3. by sample - means for figure
	# treatment * timepoint
	summ.3.a = data.inf.samp %>%
		group_by(trial, treatment, timepoint, tissue) %>%
		summarize(colonies_mean_per_g=mean(colonies_sum_per_g) ) %>%
		ungroup()

	# treatment
	summ.3.b = data.inf.samp %>%
		group_by(trial, treatment, tissue) %>%
		summarize(colonies_mean_per_g=mean(colonies_sum_per_g) ) %>%
		ungroup()
		
### 4. by sample - factor difference between timepoints
	## spread treatment means
		# make key column for spreading
		summ.4 = summ.3.a %>% mutate(timepoint_2=paste("dpi_", sprintf("%02d", timepoint), sep="") )
		
		# remove unneeded columns
		summ.4 = summ.4 %>% select(-timepoint)
		
		# spread
		summ.4 = summ.4 %>% spread(key=timepoint_2, value=colonies_mean_per_g)
	
	## calculate
	summ.4 = summ.4 %>% mutate(
		factor_diff_28v03=round(dpi_28 / dpi_03, digits=0),
		factor_diff_28v07=round(dpi_28 / dpi_07, digits=0),
		factor_diff_28v14=round(dpi_28 / dpi_14, digits=0),
		factor_diff_28v21=round(dpi_28 / dpi_21, digits=0) )
		
	## view
	summ.4 %>% arrange(tissue, treatment, trial) %>% print(n=Inf)
	
### 5. by sample - factor difference between trials 1,2 and 3,4
	## get means
	summ.5 = data.inf.samp %>%
		group_by(trial, treatment, tissue) %>%
		summarize(colonies_mean_per_g=mean(colonies_sum_per_g) ) %>%
		ungroup()
	
	## spread trial means
		# make key column for spreading
		summ.5 = summ.5 %>% mutate(trial_2=paste("trial_", trial, sep="") )
		
		# remove unneeded columns
		summ.5 = summ.5 %>% select(-trial)
		
		# spread
		summ.5 = summ.5 %>% spread(key=trial_2, value=colonies_mean_per_g)
	
	## calculate
	summ.5 = summ.5 %>% mutate(
		factor_diff_1v3=round(trial_1 / trial_3, digits=3),
		factor_diff_1v4=round(trial_1 / trial_4, digits=3),
		factor_diff_2v3=round(trial_2 / trial_3, digits=3),
		factor_diff_2v4=round(trial_2 / trial_4, digits=3) )

	## view
	summ.5 %>% arrange(tissue, treatment) %>% print(n=Inf)

### 6. by sample - factor difference between treatments (for numerical comparison)
	## get means
	summ.6 = data.inf.samp %>%
		group_by(trial, treatment, tissue) %>%
		summarize(colonies_mean_per_g=mean(colonies_sum_per_g) ) %>%
		ungroup()

	## spread by treatment		
		# make key column for spreading
		summ.6 = summ.6 %>% mutate(treatment_2=paste("trt_", treatment, sep="") )
		
		# remove unneeded columns
		summ.6 = summ.6 %>% select(-treatment)
		
		# spread
		summ.6 = summ.6 %>% spread(key=treatment_2, value=colonies_mean_per_g)
	
	## calculate
	summ.6 = summ.6 %>% mutate(
		factor_diff_crown_v_control=round(trt_crown / trt_control, digits=3),
		factor_diff_root_v_control =round(trt_root  / trt_control, digits=0),
		factor_diff_crown_v_root=round(trt_crown / trt_root, digits=3) )
		
	## view
	summ.6 %>% arrange(tissue, trial) %>% print(n=Inf)
	

###############################
#### Visualize - by sample ####
###############################
# note: counts adjusted in ggplot statement for plotting on log scale

### prepare annotations - interaction letters
	## import
	mean.letter = read_csv(file="./4_output_external/plating-roots_F_step-1-1_sas-output_mean-let - plating.csv", col_names=T, na=c("NA",""))

	## change column type
	mean.letter = mean.letter %>% mutate(trial=as.character(trial), timepoint=as.integer(timepoint) )
	
	## gather
		# gather
		mean.letter = mean.letter %>% gather(key="tissue", value="mean_letter", -trial, -treatment, -timepoint) 
		
		# remove text
		mean.letter = mean.letter %>% mutate(tissue=str_replace(tissue, "letter_", "") )
		
	## join to summary
	summ.3.a = summ.3.a %>% left_join(mean.letter, by=c(c("trial","treatment","timepoint","tissue")) )

	## make label
	summ.3.a = summ.3.a %>% mutate(label_text=case_when(
			(!is.na(mean_letter) ) ~ paste(round(colonies_mean_per_g, digits=0), "\n", mean_letter, sep=""),
			( is.na(mean_letter) ) ~ as.character(round(colonies_mean_per_g, digits=0)) ) )

### prepare annotations - treatment letters
	## make tibble
	mean.letter.trt = tibble(
		trial		= c("3","3","3"), 
		treatment   = c("control","crown","root"),
		tissue		= c("crown","crown","crown"),
		mean_letter = c("b","a","b") )
		
	## join to summary
	summ.3.b = summ.3.b %>% left_join(mean.letter.trt, by=c(c("trial","treatment","tissue")) )

### prepare
	## change order of timepoint
	data.inf.samp = data.inf.samp %>% mutate(timepoint=fct_relevel(as.character(timepoint), c("3","7","14","21","28")) )
	summ.3.a = summ.3.a %>% mutate(timepoint=fct_relevel(as.character(timepoint), c("3","7","14","21","28")) )
	
	## facet labels
	facet.lab.trial = c('1'="Trial 1", '2'="Trial 2", '3'="Trial 3", '4'="Trial 4", 'test'="Test")
	facet.lab.trt = c('control'="Control", 'crown'="Crown Inoculation", 'root'="Root Inoculation")
	
	## set values of variables
	value.adj = 0.1
	base.size = 9

### roots
	plot.roots = data.inf.samp %>% filter(tissue == "roots") %>% {
	ggplot(., aes(x=timepoint) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=2, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data={summ.3.a %>% filter(tissue == "roots") }, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
#		geom_text(data={summ.3.a %>% filter(tissue == "roots") }, aes(y=(colonies_mean_per_g+value.adj)*25, label=label_text), size=2.75, color="red", lineheight=1 ) +
		geom_text(data={summ.3.a %>% filter(tissue == "roots") }, aes(y=(colonies_mean_per_g+value.adj)*25, label=mean_letter), size=3, color="red" ) +
		facet_grid(cols=vars(trial), rows=vars(treatment), labeller=labeller(trial=facet.lab.trial, treatment=facet.lab.trt) ) +
		scale_y_continuous(trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title=element_text(size=base.size+2), axis.text=element_text(size=base.size) ) +
		theme(strip.text=element_text(size=base.size+1) ) +
		labs(x="Sample timepoint (days post inoculation)", y="Amount of M. phaseolina (colonies/g)")
	}
	ggplot2::ggsave(file="./4_output/plating_01_sample_roots.png", device="png", plot=plot.roots, width=7, height=6.5, units="in", dpi=600)
	ggplot2::ggsave(file="./4_output/Pedroncelli strawberry Mac infection fig-1.tif", device="tiff", plot=plot.roots, width=7, height=6.5, units="in", dpi=500, compression="lzw")

### crown
	## trials 1,2,4
	plot.crown.124 = data.inf.samp %>% filter(tissue == "crown" & trial %in% c("1","2","4") ) %>% {
	ggplot(., aes(x=timepoint) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=2, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data={summ.3.a %>% filter(tissue == "crown" & trial %in% c("1","2","4") ) }, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
#		geom_text(data={summ.3.a %>% filter(tissue == "crown" & trial %in% c("1","2","4") ) }, aes(y=(colonies_mean_per_g+value.adj)*25, label=label_text), size=2.75, color="red", lineheight=1 ) +
		geom_text(data={summ.3.a %>% filter(tissue == "crown" & trial %in% c("1","2","4") ) }, aes(y=(colonies_mean_per_g+value.adj)*25, label=mean_letter), size=3, color="red" ) +
		facet_grid(cols=vars(trial), rows=vars(treatment), labeller=labeller(trial=facet.lab.trial, treatment=facet.lab.trt) ) +
		scale_y_continuous(trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title=element_text(size=base.size+2), axis.text=element_text(size=base.size) ) +
		theme(strip.text=element_text(size=base.size+1) ) +
		labs(x="Sample timepoint (days post inoculation)", y="Amount of M. phaseolina (colonies/g)")
	}
	
	### trial 3
	plot.crown.3 = data.inf.samp %>% filter(tissue == "crown" & trial == "3") %>% {
	ggplot(., aes(x=treatment) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=2, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data={summ.3.b %>% filter(tissue == "crown" & trial == "3") }, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
		geom_text(data={summ.3.b %>% filter(tissue == "crown" & trial == "3") }, aes(y=(colonies_mean_per_g+value.adj)*43, label=mean_letter), size=3, color="red" ) +
		facet_grid(cols=vars(trial), labeller=labeller(trial=facet.lab.trial) ) +
		scale_y_continuous(trans="log10", limits=c(0.1,36307.81), expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title.x=element_text(size=base.size+2), axis.title.y=element_blank(), axis.text.x=element_text(size=base.size), axis.text.y=element_blank() ) +
		theme(strip.text=element_text(size=base.size+1) ) +
		labs(x="Inoculation", y="Amount of M. phaseolina (colonies/g)")
	}
	
	plot.crown.r = wrap_plots(plot_spacer(), plot.crown.3, plot_spacer(), ncol=1, heights=c(0.9,1.1,0.825) )
	
	plot.crown = wrap_plots(plot.crown.124, plot.crown.r, nrow=1, widths=c(3.1,1) )
	
	ggplot2::ggsave(file="./4_output/plating_01_sample_crown.png", device="png", plot=plot.crown, width=7, height=6.5, units="in", dpi=600)
	ggplot2::ggsave(file="./4_output/Pedroncelli strawberry Mac infection fig-2.tif", device="tiff", plot=plot.crown, width=7, height=6.5, units="in", dpi=500, compression="lzw")


##############################
#### Wash Test - Organize ####
##############################

### upload
	in.test = read_csv(file="./2_data/Inoculum Washing Test.csv", col_names=T, na=c(""))

### organize, clean up df
	## remove unneeded column
	in.test = in.test %>% select(-`Colonies/g`)

	## change column names
	in.test = in.test %>% rename(treatment_old=Treatment, plant=Rep, tissue_old=Tissue, g_tissue_plate=Weight, plate=Plate, n_colonies=Colonies)

	## reduce treatment names, lower case tissue
		# make new column
		in.test = in.test %>% mutate(
			treatment=case_when(
				(treatment_old == "Control") ~ "control",
				(treatment_old == "Crown Inoc") ~ "crown",
				(treatment_old == "Root Inoc") ~ "root"),
			tissue=case_when(
				(tissue_old == "Crown") ~ "crown",
				(tissue_old == "Roots") ~ "roots") )
		
		# remove old column
		in.test = in.test %>% select(-treatment_old, -tissue_old)

	## change column type
	in.test = in.test %>% mutate(n_colonies=as.integer(n_colonies) )
	
	## add columns to align with main df
	in.test = in.test %>% mutate(trial="test", timepoint=as.character("0") )
	
### organize - per plate	
	## calculate colonies per g
	data.test.plate = in.test %>% mutate(colonies_per_g=round(n_colonies / g_tissue_plate, digits=0) )

	## reorder columns
	data.test.plate = data.test.plate %>% select(trial, treatment, timepoint, plant, tissue, plate, g_tissue_plate, n_colonies, colonies_per_g)

### organize - per sample
	## checks
		# for missing plates
		in.test %>% group_by(trial, treatment, timepoint, plant, tissue) %>% summarize(ct=n(), ct_plate=n_distinct(plate) ) %>% filter(ct != 7 | is.na(ct) | ct_plate != 7 | is.na(ct_plate) )

		# one g_tissue_total per sample
		in.test %>% group_by(trial, treatment, timepoint, plant, tissue) %>% summarize(n_g_tissue=n_distinct(g_tissue_plate) ) %>% filter(n_g_tissue > 1)
	
	## sum colonies for each sample
	data.test.samp = in.test %>% 
		group_by(trial, treatment, timepoint, plant, tissue) %>%
		summarize(n_colonies_sum=sum(n_colonies), g_tissue_total=sum(g_tissue_plate) ) %>%
		ungroup()
		
	## calculate colonies per g
	data.test.samp = data.test.samp %>% mutate(colonies_sum_per_g=round(n_colonies_sum / g_tissue_total, digits=0) )
	
	## order columns
	data.test.samp = data.test.samp %>% select(trial, treatment, timepoint, plant, tissue, g_tissue_total, n_colonies_sum, colonies_sum_per_g)	
		
	## sort
	data.test.samp = data.test.samp %>% arrange(trial, treatment, timepoint, plant, tissue)
	
### export
	write_csv(data.test.plate, file="./2_data_curated/sb-mac-inf_data_wash-test_plate.csv", col_names=T, na=".", append=F)	
	write_csv(data.test.samp, file="./2_data_curated/sb-mac-inf_data_wash-test_sample.csv", col_names=T, na=".", append=F)	


###############################
#### Wash Test - Summarize ####
###############################

### 3. by sample - means for figure
	# treatment
	summ.test.3 = data.test.samp %>%
		group_by(trial, treatment, tissue, timepoint) %>%
		summarize(colonies_mean_per_g=mean(colonies_sum_per_g) ) %>%
		ungroup()

### 5. by sample - factor difference between trials 1,2 and 3,4 and wash
	# treatment
	summ.5 %>%
		left_join(summ.test.3) %>%
		rename(test=colonies_mean_per_g) %>%
		select(
			treatment, tissue, trial_1, trial_2, trial_3, trial_4, test, 
			factor_diff_1v3, factor_diff_1v4, factor_diff_2v3, factor_diff_2v4, 
			-trial, -timepoint) %>%
		arrange(tissue, treatment)


###############################
#### Wash Test - Visualize ####
###############################
# note: counts adjusted in ggplot statement for plotting on log scale

### prepare	
	## set values of variables
	value.adj = 0.1
	base.size = 9

### plot
	plot.test = ggplot(data.test.samp, aes(x=treatment) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=2, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data=summ.test.3, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
		facet_grid(rows=vars(tissue) ) +
		scale_y_continuous(trans="log10", breaks=c(1,10,100), expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title=element_text(size=base.size+2), axis.text=element_text(size=base.size) ) +
		theme(strip.text=element_text(size=base.size+1) ) +
		labs(x="Inoculation Treatment", y="Amount of M. phaseolina (colonies/g)")

	ggplot2::ggsave(file="./4_output/plating_01_wash-test_sample.png", device="png", plot=plot.test, width=3.25, height=3.5, units="in", dpi=600)
	ggplot2::ggsave(file="./4_output/Pedroncelli strawberry Mac infection fig-3.tif", device="tiff", plot=plot.test, width=3.25, height=3.5, units="in", dpi=500, compression="lzw")


##############################
#### Combined - Visualize ####
##############################

### prepare
	## change to factor
	data.test.samp = data.test.samp %>% mutate(timepoint=as_factor(timepoint) )
	
	## get range to set limits
	data.inf.samp %>% group_by(tissue) %>% summarize(colonies_min=min(colonies_sum_per_g), colonies_max=max(colonies_sum_per_g) )

#### combined
#	## root - test
#	p.c.1.1 = data.test.samp %>% filter(tissue == "roots") %>% {
#	ggplot(., aes(x=trial) ) +
#		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=2, position=position_jitter(w=0.05, h=0.00) ) +
#		geom_crossbar(
#			data={summ.test.3 %>% filter(tissue == "roots" ) }, 
#			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
#			width=0.3, fatten=2, color="black") +
#		facet_grid(cols=vars(tissue), rows=vars(treatment), labeller=labeller(treatment=facet.lab.trt) ) +
#		scale_y_continuous(limits=c(0.1, 8768), trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
#		theme_bw() +
#		theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(size=base.size) ) +
#		theme(strip.text.y=element_blank(), strip.text.x=element_text(size=base.size+1) ) +
#		labs(x="Trial")
#	}
#
#	## root - experiment
#	p.c.1.2 = data.inf.samp %>% filter(tissue == "roots" & timepoint %in% c("3","28") ) %>% {
#	ggplot(., aes(x=trial) ) +
#		geom_point(
#			aes(y=colonies_sum_per_g + value.adj, shape=timepoint, color=timepoint, group=timepoint), 
#			size=2, position=position_jitterdodge(jitter.height=0.00) ) +
#		geom_crossbar(
#			data={summ.3.a %>% filter(tissue == "roots" & timepoint %in% c("3","28") ) }, 
#			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj, 
#				color=timepoint, group=timepoint), 
#			position=position_dodge(width=0.75), width=0.3, fatten=2) +
#		facet_grid(cols=vars(tissue), rows=vars(treatment), labeller=labeller(treatment=facet.lab.trt) ) +
#		scale_y_continuous(limits=c(0.1, 8768), trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
#		scale_shape_manual(values=c(2,0) ) +
#		theme_bw() +
#		theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(size=base.size) ) +
#		theme(strip.text=element_text(size=base.size+1) ) +
#		guides(shape=guide_legend(direction="horizontal"), color=guide_legend(direction="horizontal") ) +
#		labs(shape="Timepoint", color="Timepoint")
#	}
#
#	## crown - test
#	p.c.2.1 = data.test.samp %>% filter(tissue == "crown") %>% {
#	ggplot(., aes(x=trial) ) +
#		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=2, position=position_jitter(w=0.05, h=0.00) ) +
#		geom_crossbar(
#			data={summ.test.3 %>% filter(tissue == "crown" ) }, 
#			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
#			width=0.3, fatten=2, color="black") +
#		facet_grid(cols=vars(tissue), rows=vars(treatment), labeller=labeller(treatment=facet.lab.trt) ) +
#		scale_y_continuous(limits=c(0.1, 8768), trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
#		theme_bw() +
#		theme(axis.title.y=element_text(size=base.size+2), axis.title.x=element_blank(), axis.text=element_text(size=base.size) ) +
#		theme(strip.text.y=element_blank(), strip.text.x=element_text(size=base.size+1) ) +
#		labs(y="Amount of M. phaseolina (colonies/g)")
#	}
#
#	## crown - experiment
#	p.c.2.2 = data.inf.samp %>% filter(tissue == "crown" & timepoint %in% c("3","28") ) %>% {
#	ggplot(., aes(x=trial) ) +
#		geom_point(
#			aes(y=colonies_sum_per_g + value.adj, shape=timepoint, color=timepoint, group=timepoint), 
#			size=2, position=position_jitterdodge(jitter.height=0.00) ) +
#		geom_crossbar(
#			data={summ.3.a %>% filter(tissue == "crown" & timepoint %in% c("3","28") ) }, 
#			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj, 
#				color=timepoint, group=timepoint), 
#			position=position_dodge(width=0.75), width=0.3, fatten=2) +
#		facet_grid(cols=vars(tissue), rows=vars(treatment), labeller=labeller(treatment=facet.lab.trt) ) +
#		scale_y_continuous(limits=c(0.1, 8768), trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
#		scale_shape_manual(values=c(2,0) ) +
#		theme_bw() +
#		theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(size=base.size) ) +
#		theme(strip.text.y=element_blank(), strip.text.x=element_text(size=base.size+1) ) +
#		guides(shape=guide_legend(direction="horizontal"), color=guide_legend(direction="horizontal") ) +
#		labs(shape="Timepoint", color="Timepoint")
#	}
#
#
#	p.c = wrap_plots(p.c.1.1, p.c.1.2, p.c.2.1, p.c.2.2, nrow=1, widths=c(1,4,1,4), guides="collect" ) +
#		plot_annotation(theme=theme(legend.position="bottom"))		
#	
#	ggplot2::ggsave(file="./4_output/plating_01_samp_experiment-wash-combined.png", device="png", plot=p.c, width=7, height=6.5, units="in", dpi=600)

### combined - all dpi
	## facet labels
	facet.lab.tissue.2 = c('crown'="Crown", 'roots'="Roots")
	facet.lab.trial.2 = c('1'="Trial 1", '2'="Trial 2", '3'="Trial 3", '4'="Trial 4", 'test'="Wash\nTest")
	facet.lab.trt.2 = c('control'="Control", 'crown'="Crown\nInoculation", 'root'="Root\nInoculation")

	## root - test
	p.c2.1.1 = data.test.samp %>% filter(tissue == "roots") %>% {
	ggplot(., aes(x=timepoint) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=1.8, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data={summ.test.3 %>% filter(tissue == "roots" ) }, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
		facet_grid(cols=vars(trial), rows=vars(treatment), labeller=labeller(trial=facet.lab.trial.2) ) +
		scale_y_continuous(limits=c(0.1, 8768), trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title.y=element_text(size=base.size+2), axis.title.x=element_blank(), axis.text=element_text(size=base.size) ) +
		theme(strip.text.y=element_blank(), strip.text.x=element_text(size=base.size+1) ) +
		labs(y="Amount of M. phaseolina (colonies/g)")
	}

	## root - experiment
	p.c2.1.2 = data.inf.samp %>% filter(tissue == "roots") %>% {
	ggplot(., aes(x=timepoint) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=1.8, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data={summ.3.a %>% filter(tissue == "roots") }, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
		facet_grid(cols=vars(trial), rows=vars(tissue, treatment), labeller=labeller(trial=facet.lab.trial.2, tissue=facet.lab.tissue.2, treatment=facet.lab.trt.2) ) +
		scale_y_continuous(trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(size=base.size) ) +
		theme(strip.text.y=element_text(size=base.size+1, margin=margin(l=1.375, r=1.375) ), strip.text.x=element_text(size=base.size+1) ) +
		labs(x="Sample timepoint (days post inoculation)", y="Amount of M. phaseolina (colonies/g)")
	}

	## crown - test
	p.c2.2.1 = data.test.samp %>% filter(tissue == "crown") %>% {
	ggplot(., aes(x=timepoint) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=1.8, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data={summ.test.3 %>% filter(tissue == "crown" ) }, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
		facet_grid(cols=vars(trial), rows=vars(treatment), labeller=labeller(trial=facet.lab.trial.2) ) +
		scale_y_continuous(limits=c(0.1, 8768), trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title.y=element_text(size=base.size+2), axis.title.x=element_blank(), axis.text=element_text(size=base.size) ) +
		theme(strip.text=element_blank() ) +
		labs(y="Amount of M. phaseolina (colonies/g)")
	}

	## crown - experiment
	p.c2.2.2 = data.inf.samp %>% filter(tissue == "crown") %>% {
	ggplot(., aes(x=timepoint) ) +
		geom_point(aes(y=colonies_sum_per_g + value.adj), shape=1, size=1.8, position=position_jitter(w=0.05, h=0.00) ) +
		geom_crossbar(
			data={summ.3.a %>% filter(tissue == "crown") }, 
			aes(y=colonies_mean_per_g + value.adj, ymin=colonies_mean_per_g + value.adj, ymax=colonies_mean_per_g + value.adj), 
			width=0.3, fatten=2, color="red") +
		facet_grid(cols=vars(trial), rows=vars(tissue, treatment), labeller=labeller(trial=facet.lab.trial.2, tissue=facet.lab.tissue.2, treatment=facet.lab.trt.2) ) +
		scale_y_continuous(trans="log10", expand=expansion(mult=c(0.05,0.075)) ) +
		theme_bw() +
		theme(axis.title.y=element_blank(), axis.title.x=element_text(size=base.size+2), axis.text.y=element_blank(), axis.text.x=element_text(size=base.size) ) +
		theme(strip.text.y=element_text(size=base.size+1, margin=margin(l=1.375, r=1.375) ), strip.text.x=element_blank() ) +
		labs(x="Sample timepoint (days post inoculation)", y="Amount of M. phaseolina (colonies/g)")
	}
	
	p.c2.1 = wrap_plots(p.c2.1.1, p.c2.1.2, nrow=1, widths=c(1.25,15) )
	
	p.c2.2 = wrap_plots(p.c2.2.1, p.c2.2.2, nrow=1, widths=c(1.25,15) )

	p.c2 = wrap_plots(p.c2.1, p.c2.2, ncol=1)
	
	ggplot2::ggsave(file="./4_output/plating_01_combined-all_experiment-wash_roots.png", device="png", plot=p.c2, width=7, height=8, units="in", dpi=600)

