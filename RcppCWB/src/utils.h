extern Corpus *corpus;
int cleanup(int error_code);
void compress_reversed_index(Attribute *attr, char *output_fn, char *corpus_id, int debug);
void decompress_check_reversed_index(Attribute *attr, char *output_fn, char *corpus_id, int debug);
int decode_check_huff(Attribute *attr, char* corpus_id, char *fname);
int compute_code_lengths(Attribute *attr, HCD *hc, char *fname);
int do_attribute(Attribute *attr, ComponentID cid, int validate);

