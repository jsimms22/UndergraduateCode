#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc_relaxed.h>

#include "packingDNAseq.h"
#include "kmer_hash_shared.h"

#define MARGIN_MULTIPLIER fmin((float)(THREADS+1)/2, 3)

int main(int argc, char *argv[]){

    /** Declarations **/
    double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

    /** Read input **/
    upc_barrier;
    inputTime -= gettime();

    //Segment the work
    int64_t n_total_kmers = getNumKmersInUFX(argv[1]);
    int64_t n_kmers_to_process_ideal = n_total_kmers / THREADS + 1;
    int64_t start_kmer = MYTHREAD * n_kmers_to_process_ideal;
    int64_t end_kmer = (MYTHREAD+1) * n_kmers_to_process_ideal;
    if (MYTHREAD == THREADS-1) end_kmer = n_total_kmers;
    int64_t n_kmers_to_process = end_kmer - start_kmer;
    int64_t char_start_position = start_kmer * LINE_SIZE;
    int64_t chars_to_read = n_kmers_to_process * LINE_SIZE;
    unsigned char * buffer = (unsigned char *)malloc(chars_to_read * sizeof(unsigned char));
    unsigned char * buffer_cpy = (unsigned char *)malloc(n_kmers_to_process_ideal * LINE_SIZE * sizeof(unsigned char));

    printf("Process %d: Reading and creating graph for K-mers %lld - %lld\n", MYTHREAD, start_kmer, end_kmer);

    //Read the appropriate portion of the data
    FILE *input_file = fopen(argv[1], "r");
    if (fseek(input_file, char_start_position, SEEK_SET) != 0) {
	printf("Error Seeking...");
	exit(0);
    }
    if (fread(buffer, sizeof(unsigned char), chars_to_read, input_file) != chars_to_read) {
	printf("Error reading...");
	exit(0);
    }
    fclose(input_file);

    upc_barrier;
    inputTime += gettime();

    /** Graph construction **/
    constrTime -= gettime();

    //Allocate memory for start kmer offsets
    int64_t *start_kmer_offsets = (int64_t *)malloc(sizeof(int64_t) * n_kmers_to_process);
    int64_t n_start_kmers = 0;

    shared shared_hash_table_t *hashtable_global;
    shared shared_memory_heap_t *memory_heap_global;
    hashtable_global = create_shared_hash_table(THREADS, MYTHREAD, n_kmers_to_process_ideal * MARGIN_MULTIPLIER);
    memory_heap_global = create_shared_memory_heap(THREADS, MYTHREAD, n_kmers_to_process_ideal * MARGIN_MULTIPLIER);
    shared_hash_table_t *private_hashtable = (shared_hash_table_t *)&hashtable_global[MYTHREAD];
    shared_memory_heap_t *private_memory_heap = (shared_memory_heap_t *)&memory_heap_global[MYTHREAD];

    //Allocate local shared memory for k-mers to be transferred to different processes
    int nints = THREADS;
    shared int64_t *process_kmer_list_offsets_global = (shared int64_t *)upc_all_alloc(THREADS, nints * sizeof(shared int64_t));
    int64_t *process_kmer_list_offsets = (int64_t *)&process_kmer_list_offsets_global[MYTHREAD];
    memset(process_kmer_list_offsets, 0, sizeof(int64_t) * THREADS);

    //We expect 1 / THREADS portion of the kmers to go to process i, multipler to allow for a margin of error
    int64_t max_kmers_to_transfer_to_single_process = (n_kmers_to_process_ideal / THREADS) * MARGIN_MULTIPLIER;
    int64_t nchars = max_kmers_to_transfer_to_single_process * LINE_SIZE * THREADS;
    shared char *kmers_to_transfer_global = (shared char *)upc_all_alloc(THREADS, nchars * sizeof(shared char));
    char *kmers_to_transfer = (char *)&kmers_to_transfer_global[MYTHREAD];

    //Read and hash kmers to individual processes
    for (int i = 0; i < n_kmers_to_process*LINE_SIZE; i += LINE_SIZE) {
	int process_owner = hashkmer(THREADS, &buffer[i]);
	//This kmer belongs to self
	if (process_owner == MYTHREAD) {
	    add_kmer_shared(private_hashtable, private_memory_heap, &buffer[i], buffer[i+KMER_LENGTH+1], buffer[i+KMER_LENGTH+2]);
	}
	//This kmer is to be transfered
	else {
	    memcpy(&kmers_to_transfer[process_owner * max_kmers_to_transfer_to_single_process * LINE_SIZE + process_kmer_list_offsets[process_owner]],
		   &buffer[i],
		   sizeof(char) * LINE_SIZE);
	    process_kmer_list_offsets[process_owner] += LINE_SIZE;
	}

	//This kmer is a start kmer
	if (buffer[i+KMER_LENGTH+1] == 'F') {
	    start_kmer_offsets[n_start_kmers++] = i;
	}
    }

    //This following section could possibly be optimized by using notify / waits (notice it re-uses buffer)
    upc_barrier;

    //Read kmers that should be transferred over to self
    for (int i = 0; i < THREADS; i++) {
	if (i != MYTHREAD) {
	    int64_t n_kmers_char_to_transfer_to_self = ((shared [0] int64_t *)&process_kmer_list_offsets_global[i] + MYTHREAD)[0];
	    upc_memget(buffer_cpy,
		       ((shared [0] char *)&kmers_to_transfer_global[i]) + MYTHREAD * max_kmers_to_transfer_to_single_process * LINE_SIZE,
		       sizeof(char) * n_kmers_char_to_transfer_to_self);

	    for (int j = 0; j < n_kmers_char_to_transfer_to_self; j += LINE_SIZE) {
		add_kmer_shared(private_hashtable, private_memory_heap, &buffer_cpy[j], buffer_cpy[j+KMER_LENGTH+1], buffer_cpy[j+KMER_LENGTH+2]);
	    }
	}
    }   

    upc_barrier;

    //Free memory for more shared space... (remove if too slow or something)
    upc_all_free(kmers_to_transfer_global);
    upc_all_free(process_kmer_list_offsets_global);

    constrTime += gettime();

    /** Graph traversal **/
    traversalTime -= gettime();

    char output_file_name[15];
    memset(output_file_name, 0, sizeof(char) * 15);
    sprintf(output_file_name, "pgen_%d.out", MYTHREAD);
    FILE *output_file = fopen(output_file_name, "w");

    char contig_seq[MAXIMUM_CONTIG_SIZE];
    memset(contig_seq, 0, sizeof(char) * MAXIMUM_CONTIG_SIZE);

    //Allocate memory for subcontigs
    for (int i = 0; i < n_start_kmers; i++) {

	//Get the kmer string
	char right_ext;
	int64_t kmer_string_offset = start_kmer_offsets[i], kmer_track = 0;
	memcpy(contig_seq, &buffer[kmer_string_offset], sizeof(char) * KMER_LENGTH);

	//Get the kmer object from distributed hash
	int process_owner_seed = hashkmer(THREADS, &buffer[kmer_string_offset]);
	shared shared_kmer_t *cur_kmer;
	if (process_owner_seed == MYTHREAD)
	    cur_kmer = lookup_kmer_shared_local(private_hashtable, &buffer[kmer_string_offset]);
	else
	    cur_kmer = lookup_kmer_shared_remote(&hashtable_global[process_owner_seed], &buffer[kmer_string_offset]);
	right_ext = cur_kmer->r_ext;

	//Extend right
	while (right_ext != 'F') {
	    kmer_track++;
	    contig_seq[kmer_track + KMER_LENGTH - 1] = right_ext;
	    int process_owner = hashkmer(THREADS, &contig_seq[kmer_track]);
	    if (process_owner == MYTHREAD)
		cur_kmer = lookup_kmer_shared_local(private_hashtable, &contig_seq[kmer_track]);
	    else
		cur_kmer = lookup_kmer_shared_remote(&hashtable_global[process_owner], &contig_seq[kmer_track]);

	    right_ext = cur_kmer->r_ext;
	}

	//Ouput contig
	fprintf(output_file, "%s\n", contig_seq);

	//Reset
	memset(contig_seq, 0, sizeof(char) * MAXIMUM_CONTIG_SIZE);
    }

    fclose(output_file);

    upc_barrier;
    traversalTime += gettime();

    /** Print timing and output info **/
    /***** DO NOT CHANGE THIS PART ****/
    if(MYTHREAD==0){
	printf("%s: Input set: %s\n", argv[0], argv[1]);
	printf("Number of UPC threads: %d\n", THREADS);
	printf("Input reading time: %f seconds\n", inputTime);
	printf("Graph construction time: %f seconds\n", constrTime);
	printf("Graph traversal time: %f seconds\n", traversalTime);
    }
    return 0;
}
