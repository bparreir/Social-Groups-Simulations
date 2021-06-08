#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "sgfunctions.h"
#include "sampling.h"
#include "mix.h"
#include "sort.h"
#include "social_groups_functions19.h"

int main(int ac, char **av) {

	/*strutures variables*/
	int s, r, t, i, j, k;
	family *population;
	indmale *aux_mlist;
	indfemale *aux_flist;
	int counter_del;
	int **nsampled; //nb of inds of each class actually sampled in a given pop
	int **idx_inds_sampled;
	int total_nsam; //total number of inds we want to sample per group
	int aux_ntotal;
	int nreps_suc;
	int **msats;
	int *nb_alleles;
	int **alleles;
	int total_nsampled;
	int **allfreq;
	double **sum_het_expect;
	double **sum_het_obs;
	double **sum_sqrthet;
	int **nrepeats_stats; //for each sampled time nb of repeats were stats were calculated
	int counter;
	double **aux_het_expct;
	double **aux_het_obs;
	double aux_het_pop;

	int nl_init; // number of lines that the file of the allelic freqs dist has
	double *aux1_h; //used to print the mean expect het of all loci per pop
	double *aux2_h; //used to print the mean observed het of all loci per pop
	int *ntotal_pop;
	double aux_gtime; //used to calculate the generation time. Sum ages of males/ females parents
	int aux_mgtime, aux_fgtime; //aux var that counts the total nb of

	/*settings variables*/
	int nfamilies;
	double gpinit_rate;
	int nfemales;
	int nmales;
	int ndomf;
	int ndomm;
	int lifet_f;
	int lifet_m;
	int max_lifet;
	int reprod_agem;
	int reprod_agef;
	int wng_age;
	int birth_interval;
	double juv_drate;
	int mean_noffspg;
	int dominance_coeff;
	int nloci;
	int ttotal;
	int tstats;
	int npsam;
	double *mut_rate;
	int *idx_psam;
	double aux_mrate;
	int aux_psam;
	int *nsam; //nb of samples from each class
	int nsim;
	int nrepeats; //nb of retetions for each param combination
	int aux_abc1;
	double aux_abc2;

	int *init_nall; //total number of alleles for each locus
	int aux_init;
	double aux_init2;
	int **init_alleles;
	double **init_allfreq;
	char hp [50];
	char dmg [50];
	char msts [50];
	char parents_age[50];

	int **cstruct;
	int **mstruct;
	int **fstruct;

	int aux_link;

	FILE *parmfile, *abcfile,*seedfile, *initfile, *clinksfile, *mlinksfile, *flinksfile, *initfreqf, *initallf;
	FILE *abc_results, *demogfile, *mstsfile, *afreqfile, *hetfile,  *hetfpop, *paragf, *gentimef, *allfreq_chosen;
/*---------------------------------------------------------------------------
 Get parameters
----------------------------------------------------------------------------*/
	parmfile = fopen("settings.txt","r");
	fscanf(parmfile,"ngroups %i\n",&nfamilies);
	fscanf(parmfile,"rate_groupinit %lf\n",&gpinit_rate);
	fscanf(parmfile,"nb_reprodfemales %d\n",&nfemales);
	fscanf(parmfile,"nb_reprodmales %d\n",&nmales);
	fscanf(parmfile,"nb_fdominants %d\n", &ndomf);
	fscanf(parmfile,"nb_mdominants %d\n", &ndomm);
	fscanf(parmfile,"noffsprings %d\n",&mean_noffspg);
	fscanf(parmfile,"domn_coefficient %d\n",&dominance_coeff);
	fscanf(parmfile,"nloci %d\n", &nloci);

	mut_rate=(double*)malloc(nloci*sizeof(double));
	fscanf(parmfile,"mrate %lf",&aux_mrate);


	mut_rate[0]=aux_mrate;
	for (i=1; i<nloci; i++) {
		fscanf(parmfile,"%lf ",&aux_mrate);
		mut_rate[i]=aux_mrate;

	}
	fscanf(parmfile,"\n");
	fscanf(parmfile,"time %d\n",&ttotal);

	fscanf(parmfile,"tstatistics %d\n", &tstats);

	fscanf(parmfile,"ngroups sampled %d\n", &npsam);

	idx_psam=(int*)malloc(npsam*sizeof(int));

	fscanf(parmfile,"idx sampled groups %d",&aux_psam);
	idx_psam[0]=aux_psam;

	for (i=1; i<npsam; i++) {
		fscanf(parmfile,"%d ", &aux_psam);
		idx_psam[i]=aux_psam;

	}
	fscanf(parmfile,"\n");

	nsam=(int*)malloc(8*sizeof(int));

	fscanf(parmfile,"nb sampled RM %d\n",&nsam[0]);
	fscanf(parmfile,"nb sampled NRm %d\n", &nsam[1]);
	fscanf(parmfile,"nb sampled juvM %d\n", &nsam[2]);
	fscanf(parmfile,"nb sampled offM %d\n", &nsam[3]);
	fscanf(parmfile,"nb sampled RF %d\n", &nsam[4]);
	fscanf(parmfile,"nb sampled NRf %d\n", &nsam[5]);
	fscanf(parmfile,"nb sampled juvF %d\n", &nsam[6]);
	fscanf(parmfile,"nb sampled offF %d\n", &nsam[7]);
	fscanf(parmfile,"nbsim %d\n", &nsim);


	fclose(parmfile);


	/*---------------------------------------------------
		get allelic frequencies
	----------------------------------------------------*/
	//initfile=fopen("initfreqs.txt","r");


	initallf=fopen("allelic_size_init.txt","r");
	fscanf(initallf,"%d ",&aux_init);
	nl_init=aux_init;
	fscanf(initallf,"\n");

	init_nall=(int*)malloc(nl_init*sizeof(int));
	for (i=0; i<nl_init; i++) {
		fscanf(initallf,"%i ",&aux_init);
		init_nall[i]=aux_init;
	}
	fscanf(initallf,"/n");

	init_alleles=(int**)malloc(nl_init*sizeof(int*));
	for (i=0; i<nl_init; i++) {
		init_alleles[i]=(int*)malloc(init_nall[i]*sizeof(int));
		for (j=0; j<init_nall[i]; j++) {
			fscanf(initallf,"%d ",&aux_init);
			init_alleles[i][j]=aux_init;

		}

		fscanf(initallf,"\n");
	}
	fclose(initallf);



	//printf("\nallelic freqs are:\n");
	initfreqf=fopen("allelic_freq_init.txt","r");
	init_allfreq=(double**)malloc(nl_init*sizeof(double*));
	for (i=0; i<nl_init; i++) {
		init_allfreq[i]=(double*)malloc(init_nall[i]*sizeof(double));
		for (j=0; j<init_nall[i]; j++){
			fscanf(initfreqf,"%lf ",&aux_init2);
			init_allfreq[i][j]=aux_init2;
		//	printf("%g ",init_allfreq[i][j]);
		}
		fscanf(initfreqf,"\n");
	}

	fclose(initfreqf);



	/*-----------------------------------------------------------------
		get colonization links
	------------------------------------------------------------------*/
	clinksfile = fopen("colonization_popstruct.txt","r");

	cstruct=(int**)malloc(nfamilies*sizeof(int*));
	for (i=0; i<nfamilies; i++) cstruct[i]=(int*)malloc(nfamilies*sizeof(int));

	for (i=0 ; i<nfamilies ; i++){
		for (j=0 ; j<nfamilies ; j++){
			fscanf(clinksfile, "%i ", &aux_link) ;
			cstruct[i][j] = aux_link ;
		}
		fscanf(clinksfile,"\n");
	}
	fclose(clinksfile);

	/*-----------------------------------------------------------------
		get males migration links
	------------------------------------------------------------------*/
	mlinksfile = fopen("mmig_struct.txt","r");

	mstruct=(int**)malloc(nfamilies*sizeof(int*));
	for (i=0; i<nfamilies; i++) mstruct[i]=(int*)malloc(nfamilies*sizeof(int));


	for (i=0 ; i<nfamilies ; i++){
		for (j=0; j<nfamilies; j++){
			fscanf(mlinksfile, "%i ", &aux_link);
			mstruct[i][j]=aux_link ;
		}
		fscanf(mlinksfile,"\n");
	}
	fclose(mlinksfile);

	/*-----------------------------------------------------------------
		get females migration links
	------------------------------------------------------------------*/
	flinksfile = fopen("fmig_struct.txt","r");

	fstruct=(int**)malloc(nfamilies*sizeof(int*));
	for (i=0; i<nfamilies; i++) fstruct[i]=(int*)malloc(nfamilies*sizeof(int));


	for (i=0 ; i<nfamilies ; i++){
		for (j=0 ; j<nfamilies ; j++){
			fscanf(flinksfile, "%i ", &aux_link) ;
			fstruct[i][j] = aux_link ;

		}
		fscanf(flinksfile,"\n");
	}
	fclose(flinksfile);

	/*------------------------------------------------------------------------
		get seed
	-------------------------------------------------------------------------*/
	seedfile = fopen("seed.txt","r");
	fscanf(seedfile,"%ld\n",&idum);
	fclose(seedfile);


	/*-------------------------------- Memory allocation -------------------------------------------
	--------------------------------------------------------------------------------------*/
	nsampled=(int**)malloc(npsam*sizeof(int*));
	for (i=0; i<npsam; i++) nsampled[i]=(int*)malloc(8*sizeof(int));

	total_nsam=nsam[0]+nsam[1]+nsam[2]+nsam[3]+nsam[4]+nsam[5]+nsam[6]+nsam[7];

	idx_inds_sampled=(int**)malloc(npsam*sizeof(int*)); //index of individuals of each class sampled in one pop
	for (i=0; i<npsam; i++) idx_inds_sampled[i]=(int*)malloc(total_nsam*sizeof(int));

	msats=(int**)malloc(nloci*sizeof(int*));
	for (i=0; i<nloci; i++) msats[i]=(int*)malloc((2*total_nsam*npsam)*sizeof(int));

	nb_alleles=(int*)malloc(nloci*sizeof(int));

	alleles=(int**)malloc(nloci*sizeof(int*));
	for (i=0; i<nloci; i++) alleles[i]=(int*)malloc((2*total_nsam*npsam)*sizeof(int));

	allfreq=(int**)malloc((npsam+1)*sizeof(int*));

	sum_het_expect=(double**)malloc((nloci+1)*sizeof(double*));
	for (i=0; i<(nloci+1); i++) sum_het_expect[i]=(double*)malloc(((ttotal/tstats)*(npsam+1))*sizeof(double));

	sum_sqrthet=(double**)malloc((nloci+1)*sizeof(double*));
	for (i=0; i<(nloci+1); i++) sum_sqrthet[i]=(double*)malloc(((ttotal/tstats)*(npsam+1))*sizeof(double));

	nrepeats_stats=(int**)malloc((nfamilies+1)*sizeof(int*));
	for (i=0; i<(nfamilies+1); i++) nrepeats_stats[i]=(int*)malloc((ttotal/tstats)*sizeof(int));

	aux_het_expct=(double**)malloc(nloci*sizeof(double*));
	for (i=0; i<nloci; i++) aux_het_expct[i]=(double*)malloc((npsam+1)*sizeof(double));

	aux_het_obs=(double**)malloc(nloci*sizeof(double*));
	for (i=0; i<nloci; i++) aux_het_obs[i]=(double*)malloc((npsam+1)*sizeof(double));

	sum_het_obs=(double**)malloc((nloci+1)*sizeof(double*));
	for (i=0; i<(nloci+1); i++) sum_het_obs[i]=(double*)malloc(((ttotal/tstats)*(npsam+1))*sizeof(double));

	aux1_h=(double*)malloc((npsam+1)*sizeof(double));

	aux2_h=(double*)malloc((npsam+1)*sizeof(double));

	ntotal_pop=(int*)malloc((npsam+1)*sizeof(int));
	/*---------------------------------------------------------------------------------
		initializes
	-----------------------------------------------------------------------------------*/
	abcfile = fopen("abc_settings.txt","r");
	abc_results=fopen("result_abc.txt","w");
	fprintf(abc_results,"lftm lftf lftmax rprdagem rprdagef wngage birthintv juvdrate nrepts nsuccess\n");

  afreqfile=fopen("allfreqs_result.txt","w");
	fprintf(afreqfile,"The results below are the number of alleles and the allelic size followed by allelic frequencies at each pop and at the overall pop\n");
	fprintf(afreqfile,"Sampled pops are:");
	for (i=0; i<npsam; i++) fprintf(afreqfile," %i",idx_psam[i]);
	fprintf(afreqfile,"\n");
	fprintf(afreqfile,"\nSim Rep t NbAll AllSize AllelicFreqs\n");

	hetfile=fopen("het_mean_result.txt","w");
	fprintf(hetfile,"The results below are the average heterozygosities between all repetions\nEach line is *one tstep* and *one pop*\nFirst line is nrepeats\nnrepeats: nb of runs used to estimate the mean het value\nAt each line results appear in the following order: nrepeats tstats expected_het pop#1,locus#1 expected_het pop1,locus#2\n");
	fprintf(hetfile,"Sampled pops are:");
	for (i=0; i<npsam; i++) fprintf(hetfile," %i",idx_psam[i]);
	fprintf(hetfile,"\n\n");

	gentimef=fopen("generation_time.txt","w");
	fprintf(gentimef,"Results bellow are the generations time.\nThis matrix has nrepeats by (2*tstats)\nFor each tstep the first value represents generation time calculated using male ages and second value represents generation time using female ages\n\n");
	allfreq_chosen=fopen("allfreqs_init.txt","w");
	fprintf(allfreq_chosen,"\nResults bellow are the allelic frequencies used to initialize each run\nEach %i lines are the allelic freqs",nloci);

	for (s=0; s<nsim; s++) {

		fscanf(abcfile,"%i ", &lifet_m);
		fscanf(abcfile,"%i ", &lifet_f);
		fscanf(abcfile,"%i ",&max_lifet);
		fscanf(abcfile,"%i ", &reprod_agem);
		fscanf(abcfile,"%i ", &reprod_agef);
		fscanf(abcfile,"%i ", &wng_age);
		fscanf(abcfile,"%i ", &birth_interval);
		fscanf(abcfile,"%lf ", &juv_drate);
		fscanf(abcfile,"%i ", &nrepeats);

		fprintf(abc_results,"%i %i %i %i %i %i %i %g %i",lifet_m,lifet_f,max_lifet,reprod_agem,reprod_agef,wng_age,birth_interval,juv_drate,nrepeats);

		nreps_suc=nrepeats;


		for (i=0; i<(nloci+1); i++) {
			for (j=0; j<((ttotal/tstats)*(npsam+1)); j++) {
				sum_het_expect[i][j]=0;
				sum_sqrthet[i][j]=0;
				sum_het_obs[i][j]=0;
			}
		}

		for (i=0; i<(npsam+1); i++) {
			for (j=0; j<(int)(ttotal/tstats); j++) nrepeats_stats[i][j]=0;
		}



		for (r=0; r<nrepeats; r++) {
			sprintf(hp,"het_sim_%i_repeat_%i.txt",s,r);
			hetfpop=fopen(hp,"w");
			fprintf(hetfpop,"\nsim #%i\nrepeat %i\nparams: lftm: %i; lftf: %i; lftmax: %i; rprdagem %i; rprdagef: %i; wngage: %i; birthintv: %i; juvdrate: %g; nrepts: %i\n",s,r,lifet_m,lifet_f,max_lifet,reprod_agem,reprod_agef,wng_age,birth_interval,juv_drate,nrepeats);


			fprintf(hetfpop,"\nThe results below are the heterozygosities, epceted and observed, for each run, each pop at each sampled tstep\nEach line are the het for one sampled population at each locus\nEach column is one locus, last column are the average het for all loci\n\nPop results appear in the following order:");

			sprintf(parents_age, "gentime_data_sim_%i_repeat_%i.txt",s,r);
			paragf=fopen(parents_age,"w");
			fprintf(paragf,"\nsim #%i\nrepeat %i\nparams: lftm: %i; lftf: %i; lftmax: %i; rprdagem %i; rprdagef: %i; wngage: %i; birthintv: %i; juvdrate: %g; nrepts: %i\n",s,r,lifet_m,lifet_f,max_lifet,reprod_agem,reprod_agef,wng_age,birth_interval,juv_drate,nrepeats);
			fprintf(paragf,"\nThe results bellow are the ages of the parents of each individual at the time of born.\nFor each ind the two numbers represent age of father followed by age of mother\n\nt pop nb_inds class\n");


			for (i=0; i<npsam; i++) fprintf(hetfpop," %i",idx_psam[i]);
			fprintf(hetfpop,"overall_pop\n");
			fprintf(hetfpop,"\nrepeat #%i\nnall %i\ntime nsamples HE HO HT\n" ,r, nloci);

			sprintf(dmg,"demog_sim_%i_repeat_%i.txt",s,r);
			demogfile=fopen(dmg,"w");
			fprintf(demogfile,"This file has the number of individuals that bellong to each class in each sim and each repeat\n");
			fprintf(demogfile,"params: lftm: %i; lftf: %i; lftmax: %i; rprdagem %i; rprdagef: %i; wngage: %i; birthintv: %i; juvdrate: %g; nrepts: %i\n",lifet_m,lifet_f,max_lifet,reprod_agem,reprod_agef,wng_age,birth_interval,juv_drate,nrepeats);
			fprintf(demogfile,"\nt pop nb.nrprdm nb.nrprdf nb.rprdm nb.rprdf nb.juvm nb.juvf nb.ofspgm nb.ofspgf\n");

			sprintf(msts,"msats_sim_%i_repeat_%i.txt",s,r);
			mstsfile=fopen(msts,"w");
			fprintf(mstsfile,"Microssatelites data\nEach line are all loci for *one individual*\n1st column is allele#1 locus #1, 2nd column is allele#2 locus#1, 3rd column allele#1 locus#2...\nsex.labels: 0-female, 1-male\nclass.labels: 1-Reprod adults 2-Non-Reprod adults 3-Juvenils 4-Offsp");
			fprintf(mstsfile,"Sampled pops are:");
			for (i=0; i<npsam; i++) fprintf(mstsfile," %i",idx_psam[i]);
			fprintf(mstsfile,"\n");
			fprintf(mstsfile,"\ntstep group sex class dominance genotype\n");

			population = (family*)malloc(nfamilies*sizeof(struct family)) ;

			init(population, nfamilies, gpinit_rate, cstruct, mstruct, fstruct, nloci, nl_init, init_nall,  init_allfreq, init_alleles, nfemales, nmales, ndomm, ndomf, lifet_m, lifet_f, max_lifet, reprod_agem, reprod_agef, allfreq_chosen);


			for (t=0; t<ttotal; t++) {
				printf("\n TIME IS=%i\n",t);

			if (ndomf==0 && ndomm==0)	reproduction (population, nfamilies, wng_age, dominance_coeff, mean_noffspg, juv_drate, nloci, lifet_m, lifet_f, max_lifet, mut_rate, t);
				printf("\nreproduction was finished!");

				death(population, nfamilies, reprod_agem, reprod_agef, ndomm, ndomf, nmales, nfemales, birth_interval);

				reorganization (population, nfamilies););
				migration (population, nfamilies, reprod_agem, reprod_agef, nmales, nfemales, nloci);
				colonization (population, nfamilies, reprod_agem, reprod_agef, nmales, nfemales, nloci);

				if (tstats!=0 && (t%tstats)==0) {
				  for (i=0; i<nfamilies; i++) {

				      fprintf(paragf,"%i %i ",t,i);
				      aux_mlist=population[i].males;
				      for (j=0; j<population[i].nmales; j++) {
					  fprintf(paragf,"%i %i ",aux_mlist->age_father,aux_mlist->age_mother);
					  aux_mlist=aux_mlist->next;
				      }
				      fprintf(paragf,"\n%i %i ",t,i);
				      aux_flist=population[i].females;
				      for (j=0; j<population[i].nfemales; j++) {
					    fprintf(paragf,"%i %i ",aux_flist->age_father,aux_flist->age_mother);
					    aux_flist=aux_flist->next;
				      }
				      fprintf(paragf,"\n");
				    }

					aux_gtime=0;
					aux_mgtime=0;
					aux_fgtime=0;
					for (i=0; i<nfamilies; i++) {
						aux_mlist=population[i].males;
						for (j=0; j<population[i].nmales; j++) {
							if (aux_mlist->age_father!=0) {
								aux_gtime+=(double)aux_mlist->age_father;
								aux_mgtime++;
								aux_mlist=aux_mlist->next;
							}
						}
						aux_flist=population[i].females;
						for (j=0; j<population[i].nfemales; j++) {
							if (aux_flist->age_father!=0) {
						 	 aux_gtime+=(double)aux_flist->age_father;
						 	 aux_mgtime++;
							  aux_flist=aux_flist->next;
							}
						}

					}
					aux_ntotal=0;
					fprintf(gentimef,"%g ",((double)aux_gtime/(double)aux_mgtime));

					//generation time using mother ages
					aux_gtime=0;
					for (i=0; i<nfamilies; i++) {
					    aux_mlist=population[i].males;
					    for (j=0; j<population[i].nmales; j++) {
							if (aux_mlist->age_mother!=0) {
								aux_gtime+=(double)aux_mlist->age_mother;
								aux_fgtime++;
								aux_mlist=aux_mlist->next;
							}
					    }
					    aux_flist=population[i].females;
					    for (j=0; j<population[i].nfemales; j++) {
							if (aux_flist->age_mother!=0) {
								aux_gtime+=(double)aux_flist->age_mother;
								aux_fgtime++;
								aux_flist=aux_flist->next;
							}

					    }
					}

					fprintf(gentimef,"%g ",((double)aux_gtime/(double)aux_fgtime));


					for (i=0; i<(npsam+1); i++) {
						aux1_h[i]=0;
						aux2_h[i]=0;
					}

					demgstats(population, nfamilies, demogfile, t, reprod_agem, reprod_agef, wng_age);


					inds_sampling(population, npsam, idx_psam, nsam, nsampled, idx_inds_sampled, reprod_agem, reprod_agef, wng_age);
					total_nsampled=0;
					for (i=0; i<npsam; i++) {
						ntotal_pop[i]=0;
						for (j=0; j<8; j++) {
							ntotal_pop[i]+=nsampled[i][j];
							total_nsampled+=nsampled[i][j];
						}
					}
					ntotal_pop[npsam]=total_nsampled;

					//sets to zero the average het between all loci for each pop
					for (i=0; i<nloci; i++) {
						for (j=0; j<(npsam+1); j++) {
							aux_het_expct[i][j]=0;
							aux_het_obs[i][j]=0;
						}
					}

					if (total_nsampled!=0) {

						msats_compute(population, npsam, idx_psam, nsampled, idx_inds_sampled, msats, nloci, mstsfile, t);
						overall_alleles(msats, nloci, npsam, nsampled, nb_alleles, alleles, nrepeats_stats, (t/tstats));
						for (i=0; i<nloci; i++) {

							for (j=0; j<(npsam+1); j++) allfreq[j]=(int*)malloc(nb_alleles[i]*sizeof(int));

							allelic_frequencies (msats[i], allfreq, npsam, nsampled, nb_alleles[i], alleles[i], afreqfile, r, t, s);

							het_expected (allfreq, nb_alleles[i], nsampled, npsam, sum_het_expect, sum_sqrthet, (t/tstats), i, aux_het_expct);

							for (j=0; j<(npsam+1); j++) free(allfreq[j]);

							het_observed(msats[i], npsam, nsampled, i, (t/tstats), sum_het_obs, aux_het_obs);


						} //end of each locus stats




						//calculates expected het (all loci) at each pop
						for (i=0; i<(npsam+1); i++) {

							aux_het_pop=0;
							for (j=0; j<nloci; j++) {
								aux_het_pop+=aux_het_expct[j][i];
							}
	
							aux1_h[i]=((double)aux_het_pop/(double)nloci);

							sum_het_expect[nloci][((t/tstats)*(npsam+1))+i]+=aux1_h[i];

						}


						for (i=0; i<(npsam+1); i++) {
							aux_het_pop=0;
							for (j=0; j<nloci; j++) {
								aux_het_pop+=aux_het_obs[j][i];
							}
							aux2_h[i]=(double)((double)aux_het_pop/(double)nloci);
							sum_het_obs[nloci][((t/tstats)*(npsam+1))+i]+=(double)((double)aux_het_pop/(double)nloci);
						}

						//prints the calculated hets in a text file
						for (i=0; i<(npsam+1); i++ ) {
							fprintf(hetfpop,"%i %i",(t*tstats),ntotal_pop[i]);
							for (j=0; j<nloci; j++) {
								fprintf(hetfpop," %g %g",aux_het_expct[j][i],aux_het_obs[j][i]);
							}
							fprintf(hetfpop," %g %g",aux1_h[i],aux2_h[i]);
							fprintf(hetfpop,"\n");
						}



					}
				}

				aux_ntotal=0;
				for (i=0; i<nfamilies; i++) aux_ntotal+=(population[i].nmales+population[i].nfemales);
				if (aux_ntotal==0) {
					nreps_suc--;
					break;
				}



			} //end of tsteps

			fprintf(gentimef,"\n");

	/*------------- counts the total nb of inds in the overall pop ----------------------------------------*/


	/*------------------------ release all memory for the next simulation ---------------------------------*/
			for (i=0; i<nfamilies; i++) {
				for (j=0; j<population[i].nmales; j++) {
					delmale (&population[i].males, 0);
				}
				for (j=0; j<population[i].nfemales; j++) {
					delfemale (&population[i].females, 0);
				}
				free(population[i].males);
				free(population[i].females);
				free(population[i].clinks);
				free(population[i].mmigration);
				free(population[i].fmigration);
				//free(population[i]);
			}
			free(population);
			fclose(hetfpop);
			fclose(demogfile);
			fclose(mstsfile);
			fclose(paragf);
		} //end of the repeats


/*--------------------------------------------------------------------------------------------------------------
	Average stats between all runs
----------------------------------------------------------------------------------------------------------------*/
		for (i=0; i<(nloci+1); i++) {

			counter=0;
			for (j=0; j<(ttotal/tstats); j++) {
				for (k=0; k<(npsam+1); k++) {

					if (nrepeats_stats[k][j]!=0) {

						sum_het_expect[i][counter]=(double)(sum_het_expect[i][counter]/(double)nrepeats_stats[k][j]);
						sum_het_obs[i][counter]=(double)(sum_het_obs[i][counter])/(double)(nrepeats_stats[k][j]);
					}
					counter++;
				}
			}

		}

		counter=0;
		for (i=0; i<(ttotal/tstats); i++) {
			for (j=0; j<(npsam+1); j++) {
				fprintf(hetfile,"%i %i",nrepeats_stats[j][i],(i*tstats));
				for (k=0; k<(nloci+1); k++) {
					fprintf(hetfile," %lf",sum_het_expect[k][counter]);
					fprintf(hetfile," %lf",sum_het_obs[k][counter]);
				}
				counter++;
				fprintf(hetfile,"\n");
			}
			
		}



		fscanf(abcfile,"\n");
		fprintf(abc_results," %i \n",nreps_suc);
		fprintf(afreqfile,"\n");

	} //end of the sims

	seedfile = fopen("seed.txt","w");
	fprintf(seedfile,"%ld",idum);
	fclose(seedfile);

	fclose(abcfile);
	fclose(abc_results);
	fclose(afreqfile);
	fclose(hetfile);
	fclose(gentimef);
	fclose(allfreq_chosen);

}
