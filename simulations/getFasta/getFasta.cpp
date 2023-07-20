
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <vector>

#define MAX_LINE_LENGTH 1024
#define MAX_INDIVIDUALS 1024
#define MAX_SITES 1024


#define ASSERT(expr)                                                                                              \
	do{ \
		if (!((expr))) {                                                                                              \
			fprintf(stderr, "\n\n*******\n[ERROR](%s/%s:%d) %s\n*******\n", __FILE__, __FUNCTION__, __LINE__, #expr); \
			exit(1);                                                                                                  \
		} \
	}while(0);


int main(int argc, char *argv[]) {
	char *input_filename = NULL;
	char *output_prefix = NULL;
	int seed = 42;
	int opt;

	while ((opt = getopt(argc, argv, "i:o:s:")) != -1) {
		switch (opt) {
			case 'i':
				input_filename = optarg;
				break;
			case 'o':
				output_prefix = optarg;
				break;
			case 's':
				seed = atoi(optarg);
				break;
			default:
				printf("Usage: %s -i input_file -o output_prefix -s seed_value\n", argv[0]);
				return -1;
		}
	}


	if (input_filename == NULL || output_prefix == NULL) {
		printf("Usage: %s -i input_file -o output_prefix -s seed_value\n", argv[0]);
		return -1;
	}

	FILE *file = fopen(input_filename, "r");
	if (file == NULL) {
		printf("Could not open file: %s\n", input_filename);
		return -1;
	}

	if(seed==-1){
		fprintf(stderr,"Setting seed to 42\n");
		seed=42;
	}
	fprintf(stderr,"Setting seed to %d\n",seed);
	srand(seed);

	char line[MAX_LINE_LENGTH];
	fgets(line, MAX_LINE_LENGTH, file);  // Skip header line

	int tmp_nSites=4096;
	int tmp_nInds=40;


	/* indData[ind_i][site_i] */
	char** indData= (char**)malloc(tmp_nInds*sizeof(char*));
	for(int i=0;i<tmp_nInds;++i){
		indData[i]=(char*)malloc(tmp_nSites*sizeof(char));
		for(int x=0;x<tmp_nSites;++x){
			indData[i][x]='N';
		}
	}

	int site=-1;
	int ind=-1;
	int na=-1;
	int nc=-1;
	int ng=-1;
	int nt=-1;
	
	const char bases[4] = {'A', 'C', 'G', 'T'};

	// assume sites are sorted and start from 0
	int last_site=-1; //index
	int last_ind=-1; //index

	int nInds=0;
	int nSites=0;

	while (fgets(line, MAX_LINE_LENGTH, file)) {
		sscanf(line, "%d\t%d\t%d\t%d\t%d\t%d", &site, &ind, &na, &nc, &ng, &nt);


		if(ind>last_ind){
			nInds++;
			last_ind=ind;
		}

		if(site>last_site){
			nSites++;
			last_site=site;
		}

		ASSERT(site>=last_site);

		if(tmp_nSites<nSites){
			int old_size=tmp_nSites;
			tmp_nSites+=4096;
			for(int i=0;i<tmp_nInds;++i){
				char* t;
				t=(char*)realloc(indData[i],tmp_nSites*sizeof(char));
				ASSERT(t!=NULL);
				indData[i]=t;
				for(int x=old_size;x<tmp_nSites;++x){
					indData[i][x]='N';
				}
			}
		}

		int counts[4] = {na, nc, ng, nt};

		int max =0;
		int nmax=0;
		for (int i = 0; i < 4; ++i) {

			// new equal max
			if (max>0 && counts[i] == max) {
				nmax++;
				continue;
			}

			// brand new max
			if (counts[i] > max) {
				nmax=1;
				max = counts[i];
				continue;
			}
		}

		if(nmax==0){
			continue;
		}

		char max_bases[nmax];
		int num_max_bases = 0;
		for (int i = 0; i < 4; ++i) {
			if (counts[i] == max) {
				max_bases[num_max_bases] = bases[i];
				num_max_bases++;
			}
		}
		ASSERT(num_max_bases==nmax);
		ASSERT(num_max_bases>0);

		int rand_base_index = rand() % num_max_bases;
		const char base = max_bases[rand_base_index];

		indData[ind][site] = base;

	}
	fclose(file);

	ASSERT(tmp_nInds>=nInds);
	/* for (int i = 0; i < nInds; ++i) { */
		/* char output_filename[MAX_LINE_LENGTH]; */
		/* sprintf(output_filename, "%s_ind_%d.fasta", output_prefix, i); */
		/* FILE *output_file = fopen(output_filename, "w"); */
/*  */
		/* fprintf(stderr, "\nWriting to: %s",output_filename); */
/*  */
		/* fprintf(output_file, ">%d\n",i); */
		/* for (int j = 0; j < nSites; ++j) { */
			/* fprintf(output_file, "%c",indData[i][j]); */
		/* } */
		/* fprintf(output_file, "\n",i); */
/*  */
		/* fclose(output_file); */
	/* } */

	//printto same file
	char output_filename[MAX_LINE_LENGTH];
	sprintf(output_filename, "%s.fasta", output_prefix);
	FILE *output_file = fopen(output_filename, "w");
	fprintf(stderr, "\nWriting to: %s",output_filename);


	fprintf(stderr,"\nTotal number of sites: %d\n",nSites);
	for (int i = 0; i < nInds; ++i) {
		fprintf(output_file, ">%d\n",i);
		for(int s=0;s<nSites;++s){
			fprintf(output_file, "%c",indData[i][s]);
		}
		fprintf(output_file, "\n");
	}
	fclose(output_file);


	for(int i=0;i<tmp_nInds;++i){
		free(indData[i]);
		indData[i]=NULL;
	}
	free(indData);
	indData=NULL;


	return 0;
}
