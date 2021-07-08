/* X3F_IO.C
 *
 * Library for accessing X3F Files.
 *
 * Copyright 2015 - Roland and Erik Karlsson
 * BSD-style - see doc/copyright.txt
 *
 */
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wnarrowing"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include "x3f_io.h"
#include "x3f_printf.h"

#include <cstring>
#include <cstdlib>
#include <cassert>

#if defined(_WIN32) || defined (_WIN64)
#include <windows.h>
#else
#include <iconv.h>
#endif

using std::uint8_t;
using std::uint16_t;
using std::uint32_t;
using std::size_t;
using std::int8_t;
using std::int16_t;
using std::int32_t;
//using std::iconv_t;

#include "../myfile.h" //#include <cstdio>
//using std::FILE;
//using std::fseek;
//using std::ftell;

using std::malloc;
using std::realloc;
using std::calloc;
using std::free;
//using std::alloca;


namespace x3ftools{
    
/* --------------------------------------------------------------------- */
/* Hacky external flags                                                 */
/* --------------------------------------------------------------------- */
/* extern */ int legacy_offset = 0;
/* extern */ bool_t auto_legacy_offset = 1;

/* --------------------------------------------------------------------- */
/* Huffman Decode Macros                                                 */
/* --------------------------------------------------------------------- */
constexpr int HUF_TREE_MAX_LENGTH = 27;
inline constexpr int HUF_TREE_MAX_NODES(int _leaves){ return (HUF_TREE_MAX_LENGTH+1) * _leaves;}
inline constexpr uint32_t HUF_TREE_GET_LENGTH(uint32_t _v){ return ((_v>>27)&0x1f);}
inline constexpr uint32_t HUF_TREE_GET_CODE(uint32_t _v){ return (_v&0x07ffffff);}


/* --------------------------------------------------------------------- */
/* Reading and writing - assuming little endian in the file              */
/* --------------------------------------------------------------------- */
static inline void FREE(void *P){
    free(P); 
    P = nullptr;
}

static inline void GET_PROPERTY_TABLE(x3f_property_table_t *T, int NUM, x3f_info_t *I){
    int i;
    T->size = NUM;
    T->element = static_cast<x3f_property_t *>(realloc(T->element, NUM * sizeof(T->element[0])));
    for (i=0; i<T->size; i++){
        GET4(&T->element[i].name_offset, I);
        GET4(&T->element[i].value_offset, I);
    }
}

static inline void GET_TRUE_HUFF_TABLE(x3f_true_huffman_t *T, x3f_info_t *I){
    int i;
    T->element = nullptr;
    for (i=0; ; i++){
      T->size = i + 1;						\
      T->element = static_cast<x3f_true_huffman_element_t *>
                     (realloc(T->element, (i + 1) * sizeof(T->element[0])));
      GET1(&T->element[i].code_size, I);
      GET1(&T->element[i].code, I);
      if (T->element[i].code_size==0)
          break;
    }
}

      
/* --------------------------------------------------------------------- */
/* Allocating Huffman tree help data                                   */
/* --------------------------------------------------------------------- */
static void cleanup_huffman_tree(x3f_hufftree_t *HTP){
    free(HTP->nodes);
}


static void new_huffman_tree(x3f_hufftree_t *HTP, int bits){
    int leaves = 1<<bits;
    HTP->free_node_index = 0;
    HTP->nodes = (x3f_huffnode_t *)
        calloc(1, HUF_TREE_MAX_NODES(leaves) * sizeof(x3f_huffnode_t));
}


/* --------------------------------------------------------------------- */
/* Allocating TRUE engine RAW help data                                  */
/* --------------------------------------------------------------------- */
static void cleanup_true(x3f_true_t **TRUP){
    x3f_true_t *TRU = *TRUP;
    if (TRU==nullptr) return;

    x3f_printf(DEBUG, "Cleanup TRUE data\n");
    FREE(TRU->table.element);
    FREE(TRU->plane_size.element);
    cleanup_huffman_tree(&TRU->tree);
    FREE(TRU->x3rgb16.buf);
    FREE(TRU);
    *TRUP = nullptr;
}


static x3f_true_t *new_true(x3f_true_t **TRUP){
    x3f_true_t *TRU = (x3f_true_t *)calloc(1, sizeof(x3f_true_t));
    cleanup_true(TRUP);
    TRU->table.size = 0;
    TRU->table.element = nullptr;
    TRU->plane_size.size = 0;
    TRU->plane_size.element = nullptr;
    TRU->tree.nodes = nullptr;
    TRU->x3rgb16.data = nullptr;
    TRU->x3rgb16.buf = nullptr;
    *TRUP = TRU;
    return TRU;
}


static void cleanup_quattro(x3f_quattro_t **QP){
    x3f_quattro_t *Q = *QP;
    if (Q==nullptr) return;

    x3f_printf(DEBUG, "Cleanup Quattro\n");
    FREE(Q->top16.buf);
    FREE(Q);
    *QP = nullptr;
}


static x3f_quattro_t *new_quattro(x3f_quattro_t **QP){
    x3f_quattro_t *Q = (x3f_quattro_t *)calloc(1, sizeof(x3f_quattro_t));
    int i;

    cleanup_quattro(QP);
    for (i=0; i<TRUE_PLANES; i++){
      Q->plane[i].columns = 0;
      Q->plane[i].rows = 0;
    }
    Q->unknown = 0;
    Q->top16.data = nullptr;
    Q->top16.buf = nullptr;
    *QP = Q;
    return Q;
}


/* --------------------------------------------------------------------- */
/* Allocating Huffman engine help data                                   */
/* --------------------------------------------------------------------- */
static void cleanup_huffman(x3f_huffman_t **HUFP){
    x3f_huffman_t *HUF = *HUFP;
    if (HUF==nullptr) return;

    x3f_printf(DEBUG, "Cleanup Huffman\n");
    FREE(HUF->mapping.element);
    FREE(HUF->table.element);
    cleanup_huffman_tree(&HUF->tree);
    FREE(HUF->row_offsets.element);
    FREE(HUF->rgb8.buf);
    FREE(HUF->x3rgb16.buf);
    FREE(HUF);
    *HUFP = nullptr;
}


static x3f_huffman_t *new_huffman(x3f_huffman_t **HUFP){
    x3f_huffman_t *HUF = (x3f_huffman_t *)calloc(1, sizeof(x3f_huffman_t));

    cleanup_huffman(HUFP);
    /* Set all not read data block pointers to nullptr */
    HUF->mapping.size = 0;
    HUF->mapping.element = nullptr;
    HUF->table.size = 0;
    HUF->table.element = nullptr;
    HUF->tree.nodes = nullptr;
    HUF->row_offsets.size = 0;
    HUF->row_offsets.element = nullptr;
    HUF->rgb8.data = nullptr;
    HUF->rgb8.buf = nullptr;
    HUF->x3rgb16.data = nullptr;
    HUF->x3rgb16.buf = nullptr;
    *HUFP = HUF;
   return HUF;
}


/* --------------------------------------------------------------------- */
/* Creating a new x3f structure from file                                */
/* --------------------------------------------------------------------- */
/* extern */ 
x3f_t *x3f_new_from_file(IMFILE *infile){
    x3f_t *x3f = (x3f_t *)calloc(1, sizeof(x3f_t));
    x3f_info_t *I = nullptr;
    x3f_header_t *H = nullptr;
    x3f_directory_section_t *DS = nullptr;
    size_t d;
    int i;

    I = &x3f->info;
    I->error = nullptr;
    I->input.file = infile;
    I->output.file = nullptr;
    if (infile==nullptr){
        I->error = "No infile";
        return x3f;
    }

    /* Read file header */
    H = &x3f->header;
    fseek(infile, 0, SEEK_SET);
    GET4(&H->identifier, I);

    if (H->identifier!=X3F_FOVb){
        x3f_printf(ERR, "Faulty file type\n");
        x3f_delete(x3f);
        return nullptr;
    }

    GET4(&H->version, I);
    GETN(H->unique_identifier, SIZE_UNIQUE_IDENTIFIER, I);
    
    /* TODO: the meaning of the rest of the header for version >= 4.0
     * (Quattro) is unknown */
    if (H->version<X3F_VERSION_4_0){
        GET4(&H->mark_bits, I);
        GET4(&H->columns, I);
        GET4(&H->rows, I);
        GET4(&H->rotation, I);
        if (H->version>=X3F_VERSION_2_1){
            int num_ext_data = H->version >= X3F_VERSION_3_0 ? 
                               NUM_EXT_DATA_3_0 : NUM_EXT_DATA_2_1;
            GETN(H->white_balance, SIZE_WHITE_BALANCE, I);
            if (H->version>=X3F_VERSION_2_3)
                GETN(H->color_mode, SIZE_COLOR_MODE, I);
            GETN(H->extended_types, num_ext_data, I);
            for (i=0; i<num_ext_data; i++)
                GET4F(&H->extended_data[i], I);
        }
    }

    /* Go to the beginning of the directory */
    fseek(infile, -4, SEEK_END);
    fseek(infile, x3f_get4(infile), SEEK_SET);

    /* Read the directory header */
    DS = &x3f->directory_section;
    GET4(&DS->identifier, I);
    GET4(&DS->version, I);
    GET4(&DS->num_directory_entries, I);

    if (DS->num_directory_entries>0){
        size_t size = DS->num_directory_entries * sizeof(x3f_directory_entry_t);
        DS->directory_entry = (x3f_directory_entry_t *)calloc(1, size);
    }

    /* Traverse the directory */
    for (d=0; d<DS->num_directory_entries; d++){
        x3f_directory_entry_t *DE = &DS->directory_entry[d];
        x3f_directory_entry_header_t *DEH = &DE->header;
        uint32_t save_dir_pos;

        /* Read the directory entry info */
        GET4(&DE->input.offset, I);
        GET4(&DE->input.size, I);

        DE->output.offset = 0;
        DE->output.size = 0;

        GET4(&DE->_type, I);

        /* Save current pos and go to the entry */
        save_dir_pos = ftell(infile);
        fseek(infile, DE->input.offset, SEEK_SET);

        /* Read the type independent part of the entry header */
        DEH = &DE->header;
        GET4(&DEH->identifier, I);
        GET4(&DEH->version, I);

        /* NOTE - the tests below could be made on DE->_type instead */
        if (DEH->identifier==X3F_SECp){
            x3f_property_list_t *PL = &DEH->data_subsection.property_list;

            /* Read the property part of the header */
            GET4(&PL->num_properties, I);
            GET4(&PL->character_format, I);
            GET4(&PL->reserved, I);
            GET4(&PL->total_length, I);

            /* Set all not read data block pointers to nullptr */
            PL->data = nullptr;
            PL->data_size = 0;
        }

        if (DEH->identifier == X3F_SECi) {
            x3f_image_data_t *ID = &DEH->data_subsection.image_data;

            /* Read the image part of the header */
            GET4(&ID->_type, I);
            GET4(&ID->format, I);
            ID->type_format = (ID->_type << 16) + (ID->format);
            GET4(&ID->columns, I);
            GET4(&ID->rows, I);
            GET4(&ID->row_stride, I);

            /* Set all not read data block pointers to nullptr */
            ID->huffman = nullptr;

            ID->data = nullptr;
            ID->data_size = 0;
        }

        if (DEH->identifier == X3F_SECc) {
            x3f_camf_t *CAMF = &DEH->data_subsection.camf;

            /* Read the CAMF part of the header */
            GET4(&CAMF->_type, I);
            GET4(&CAMF->tN.val0, I);
            GET4(&CAMF->tN.val1, I);
            GET4(&CAMF->tN.val2, I);
            GET4(&CAMF->tN.val3, I);

            /* Set all not read data block pointers to nullptr */
            CAMF->data = nullptr;
            CAMF->data_size = 0;

            /* Set all not allocated help pointers to nullptr */
            CAMF->table.element = nullptr;
            CAMF->table.size = 0;
            CAMF->tree.nodes = nullptr;
            CAMF->decoded_data = nullptr;
            CAMF->decoded_data_size = 0;
            CAMF->entry_table.element = nullptr;
            CAMF->entry_table.size = 0;
        }
        /* Reset the file pointer back to the directory */
        fseek(infile, save_dir_pos, SEEK_SET);
    }
    return x3f;
} // x3f_new_from_file


/* --------------------------------------------------------------------- */
/* Clean up an x3f structure                                             */
/* --------------------------------------------------------------------- */
static void free_camf_entry(camf_entry_t *entry){
    FREE(entry->property_name);
    FREE(entry->property_value);
    FREE(entry->matrix_decoded);
    FREE(entry->matrix_dim_entry);
}


/* extern */ 
x3f_return_t x3f_delete(x3f_t *x3f){
    x3f_directory_section_t *DS;
    size_t d;

    if (x3f==nullptr)
        return X3F_ARGUMENT_ERROR;
    x3f_printf(DEBUG, "X3F Delete\n");

    DS = &x3f->directory_section;

    for (d=0; d<DS->num_directory_entries; d++){
        x3f_directory_entry_t *DE = &DS->directory_entry[d];
        x3f_directory_entry_header_t *DEH = &DE->header;

        if (DEH->identifier == X3F_SECp) {
            x3f_property_list_t *PL = &DEH->data_subsection.property_list;
            size_t i;

            for (i=0; i<PL->property_table.size; i++) {
                FREE(PL->property_table.element[i].name_utf8);
                FREE(PL->property_table.element[i].value_utf8);
            }
            FREE(PL->property_table.element);
            FREE(PL->data);
        }

        if (DEH->identifier==X3F_SECi){
            x3f_image_data_t *ID = &DEH->data_subsection.image_data;
            cleanup_huffman(&ID->huffman);
            cleanup_true(&ID->tru);
            cleanup_quattro(&ID->quattro);
            FREE(ID->data);
        }

        if (DEH->identifier==X3F_SECc){
            x3f_camf_t *CAMF = &DEH->data_subsection.camf;
            size_t i;

            FREE(CAMF->data);
            FREE(CAMF->table.element);
            cleanup_huffman_tree(&CAMF->tree);
            FREE(CAMF->decoded_data);
            for (i=0; i<CAMF->entry_table.size; i++){
                free_camf_entry(&CAMF->entry_table.element[i]);
            }
            FREE(CAMF->entry_table.element);
        }
    }
    FREE(DS->directory_entry);
    FREE(x3f);
    return X3F_OK;
} // x3f_delete


/* --------------------------------------------------------------------- */
/* Getting a reference to a directory entry                              */
/* --------------------------------------------------------------------- */
/* TODO: all those only get the first instance */
static x3f_directory_entry_t *x3f_get(x3f_t *x3f, uint32_t _type, uint32_t image_type){
    x3f_directory_section_t *DS;
    size_t d;

    if (x3f==nullptr) 
        return nullptr;

    DS = &x3f->directory_section;

    for (d=0; d<DS->num_directory_entries; d++){
        x3f_directory_entry_t *DE = &DS->directory_entry[d];
        x3f_directory_entry_header_t *DEH = &DE->header;
        if (DEH->identifier==_type){
            switch (DEH->identifier){
                case X3F_SECi:
                {
                    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
                    if (ID->type_format==image_type)
                        return DE;
                }
                break;
                default:
                    return DE;
            }
        }
    }
    return nullptr;
}


/* extern */ 
x3f_directory_entry_t *x3f_get_raw(x3f_t *x3f){
    x3f_directory_entry_t *DE;
    if ((DE = x3f_get(x3f, X3F_SECi, X3F_IMAGE_RAW_HUFFMAN_X530)) != nullptr)
        return DE;
    if ((DE = x3f_get(x3f, X3F_SECi, X3F_IMAGE_RAW_HUFFMAN_10BIT)) != nullptr)
        return DE;
    if ((DE = x3f_get(x3f, X3F_SECi, X3F_IMAGE_RAW_TRUE)) != nullptr)
        return DE;
    if ((DE = x3f_get(x3f, X3F_SECi, X3F_IMAGE_RAW_MERRILL)) != nullptr)
        return DE;
    if ((DE = x3f_get(x3f, X3F_SECi, X3F_IMAGE_RAW_QUATTRO)) != nullptr)
        return DE;
    if ((DE = x3f_get(x3f, X3F_SECi, X3F_IMAGE_RAW_SDQ)) != nullptr)
        return DE;
    if ((DE = x3f_get(x3f, X3F_SECi, X3F_IMAGE_RAW_SDQH)) != nullptr)
        return DE;
    return nullptr;
}


/* extern */ 
x3f_directory_entry_t *x3f_get_thumb_plain(x3f_t *x3f){
    return x3f_get(x3f, X3F_SECi, X3F_IMAGE_THUMB_PLAIN);
}


/* extern */ 
x3f_directory_entry_t *x3f_get_thumb_huffman(x3f_t *x3f){   
    return x3f_get(x3f, X3F_SECi, X3F_IMAGE_THUMB_HUFFMAN);
}


/* extern */ 
x3f_directory_entry_t *x3f_get_thumb_jpeg(x3f_t *x3f){
    return x3f_get(x3f, X3F_SECi, X3F_IMAGE_THUMB_JPEG);
}


/* extern */ 
x3f_directory_entry_t *x3f_get_camf(x3f_t *x3f){
    return x3f_get(x3f, X3F_SECc, 0);
}


/* extern */ 
x3f_directory_entry_t *x3f_get_prop(x3f_t *x3f){
    return x3f_get(x3f, X3F_SECp, 0);
}


/* For some obscure reason, the bit numbering is weird. It is
 * generally some kind of "big endian" style - e.g. the bit 7 is the
 * first in a byte and bit 31 first in a 4 byte int. For patterns in
 * the huffman pattern table, bit 27 is the first bit and bit 26 the
 * next one. */
static inline int PATTERN_BIT_POS(int len, int bit){
    return (len - bit - 1);
}

static inline int MEMORY_BIT_POS(int bit){
    return PATTERN_BIT_POS(8, bit);
}


/* --------------------------------------------------------------------- */
/* Huffman Decode                                                        */
/* --------------------------------------------------------------------- */
/* Make the huffman tree */

#ifdef DBG_PRNT
static char *display_code(int length, uint32_t code, char *buffer){
    int i;
    for (i=0; i<length; i++){
        int pos = PATTERN_BIT_POS(length, i);
        buffer[i] = ((code>>pos)&1) == 0 ? '0' : '1';
    }
    buffer[i] = 0;
    return buffer;
}
#endif


static x3f_huffnode_t *new_node(x3f_hufftree_t *tree){
    x3f_huffnode_t *t = &tree->nodes[tree->free_node_index];
    t->branch[0] = nullptr;
    t->branch[1] = nullptr;
    t->leaf = UNDEFINED_LEAF;
    tree->free_node_index++;
    return t;
}


static void add_code_to_tree(x3f_hufftree_t *tree, int length, 
                             uint32_t code, uint32_t value){
    int i;
    x3f_huffnode_t *t = tree->nodes;

    for (i=0; i<length; i++){
        int pos = PATTERN_BIT_POS(length, i);
        int bit = (code>>pos)&1;
        x3f_huffnode_t *t_next = t->branch[bit];
        if (t_next==nullptr)
            t_next = t->branch[bit] = new_node(tree);
        t = t_next;
    }
    t->leaf = value;
}


static void populate_true_huffman_tree(x3f_hufftree_t *tree, x3f_true_huffman_t *table){
    size_t i;
    new_node(tree);

    for (i=0; i<table->size; i++){
        x3f_true_huffman_element_t *element = &table->element[i];
        uint32_t length = element->code_size;

        if (length!=0){
            /* add_code_to_tree wants the code right adjusted */
            uint32_t code = ((element->code) >> (8 - length)) & 0xff;
            uint32_t value = i;
            add_code_to_tree(tree, length, code, value);

#ifdef DBG_PRNT
            {
                char buffer[100];
                x3f_printf(DEBUG, "H %5d : %5x : %5d : %02x %08x (%08x) (%s)\n",
                       i, i, value, length, code, value, display_code(length, code, buffer));
            }
#endif
        }   
    }
}


static void populate_huffman_tree(x3f_hufftree_t *tree,
				  x3f_table32_t *table, x3f_table16_t *mapping){
    size_t i;
    new_node(tree);

    for (i=0; i<table->size; i++){
        uint32_t element = table->element[i];

        if (element!=0){
            uint32_t length = HUF_TREE_GET_LENGTH(element);
            uint32_t code = HUF_TREE_GET_CODE(element);
            uint32_t value;

            /* If we have a valid mapping table - then the value from the
             mapping table shall be used. Otherwise we use the current
             index in the table as value. */
            if (table->size == mapping->size)
                value = mapping->element[i];
            else
                value = i;
            add_code_to_tree(tree, length, code, value);

#ifdef DBG_PRNT
            {
                char buffer[100];
                x3f_printf(DEBUG, "H %5d : %5x : %5d : %02x %08x (%08x) (%s)\n",
                   i, i, value, length, code, element, display_code(length, code, buffer));
            }
#endif
        }
    }
}


#ifdef DBG_PRNT
static void print_huffman_tree(x3f_huffnode_t *t, int length, uint32_t code){
    char buf1[100];
    char buf2[100];

    x3f_printf(DEBUG, "%*s (%s,%s) %s (%s)\n",
	     length, length < 1 ? "-" : (code&1) ? "1" : "0",
	     t->branch[0]==nullptr ? "-" : "0",
	     t->branch[1]==nullptr ? "-" : "1",
	     t->leaf==UNDEFINED_LEAF ? "-" : (sprintf(buf1, "%x", t->leaf),buf1),
	     display_code(length, code, buf2));

      code = code << 1;
      if (t->branch[0]) print_huffman_tree(t->branch[0], length+1, code+0);
      if (t->branch[1]) print_huffman_tree(t->branch[1], length+1, code+1);
}
#endif


/* Help machinery for reading bits in a memory */
typedef struct bit_state_s{
      uint8_t *next_address;
      uint8_t bit_offset;
      uint8_t bits[8];
}bit_state_t;


static void set_bit_state(bit_state_t *BS, uint8_t *address){
      BS->next_address = address;
      BS->bit_offset = 8;
}


static uint8_t get_bit(bit_state_t *BS){
    if (BS->bit_offset==8){
        uint8_t _byte = *BS->next_address;
        int i;

        for (i=7; i>=0; i--){
            BS->bits[i] = _byte&1;
            _byte = _byte >> 1;
        }
        BS->next_address++;
        BS->bit_offset = 0;
    }
    return BS->bits[BS->bit_offset++];
}


/* Decode use the TRUE algorithm */
static int32_t get_true_diff(bit_state_t *BS, x3f_hufftree_t *HTP){
    int32_t diff;
    x3f_huffnode_t *node = &HTP->nodes[0];
    uint8_t bits;

    while (node->branch[0]!=nullptr || node->branch[1]!=nullptr){
        uint8_t bit = get_bit(BS);
        x3f_huffnode_t *new_node = node->branch[bit];
        node = new_node;
        if (node==nullptr){
            /* TODO: Shouldn't this be treated as a fatal error? */
            x3f_printf(ERR, "Huffman coding got unexpected bit\n");
            return 0;
        }
    }

    bits = node->leaf;
    if (bits==0)
        diff = 0;
    else {
        uint8_t first_bit = get_bit(BS);
        int i;
        diff = first_bit;
        for (i=1; i<bits; i++)
            diff = (diff << 1) + get_bit(BS);
        if (first_bit==0)
            diff -= (1<<bits) - 1;
    }
    return diff;
}


/* This code (that decodes one of the X3F color planes, really is a
   decoding of a compression algorithm suited for Bayer CFA data. In
   Bayer CFA the data is divided into 2x2 squares that represents
   (R,G1,G2,B) data. Those four positions are (in this compression)
   treated as one data stream each, where you store the differences to
   previous data in the stream. The reason for this is, of course,
   that the date is more often than not near to the next data in a
   stream that represents the same color. */
/* TODO: write more about the compression */

static void true_decode_one_color(x3f_image_data_t *ID, int color){
    x3f_true_t *TRU = ID->tru;
    x3f_quattro_t *Q = ID->quattro;
    uint32_t seed = TRU->seed[color]; /* TODO : Is this correct ? */
    size_t row;
    x3f_hufftree_t *tree = &TRU->tree;
    bit_state_t BS;
    int32_t row_start_acc[2][2];
    uint32_t rows = ID->rows;
    uint32_t cols = ID->columns;
    x3f_area16_t *area = &TRU->x3rgb16;
    uint16_t *dst = area->data + color;

    set_bit_state(&BS, TRU->plane_address[color]);

    row_start_acc[0][0] = seed;
    row_start_acc[0][1] = seed;
    row_start_acc[1][0] = seed;
    row_start_acc[1][1] = seed;

    if (ID->type_format == X3F_IMAGE_RAW_QUATTRO ||
            ID->type_format == X3F_IMAGE_RAW_SDQ ||
            ID->type_format == X3F_IMAGE_RAW_SDQH){
        rows = Q->plane[color].rows;
        cols = Q->plane[color].columns;
        if (Q->quattro_layout && color==2){
            area = &Q->top16;
            dst = area->data;
        }
        x3f_printf(DEBUG, "Quattro decode one color (%d) rows=%d cols=%d\n",
           color, rows, cols);
    } else {
        x3f_printf(DEBUG, "TRUE decode one color (%d) rows=%d cols=%d\n",
                   color, rows, cols);
    }

    assert(rows==area->rows && cols>=area->columns);

    for (row = 0; row<rows; row++){
        size_t col;
        bool_t odd_row = row&1;
        int32_t acc[2];

        for (col = 0; col<cols; col++){
            bool_t odd_col = col&1;
            int32_t diff = get_true_diff(&BS, tree);
            int32_t prev = col < 2 ?
            row_start_acc[odd_row][odd_col] :
            acc[odd_col];
            int32_t value = prev + diff;

            acc[odd_col] = value;
            if (col<2)
                row_start_acc[odd_row][odd_col] = value;

            /* Discard additional data at the right for binned Quattro plane 2 */
            if (col>=area->columns) continue;

            *dst = value;
            dst += area->channels;
        }
    }
} // true_decode_one_color


static void true_decode(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    int color;
    for (color=0; color<3; color++){
        true_decode_one_color(ID, color);
    }
}


/* Decode use the huffman tree */
static int32_t get_huffman_diff(bit_state_t *BS, x3f_hufftree_t *HTP){
    int32_t diff;
    x3f_huffnode_t *node = &HTP->nodes[0];

    while (node->branch[0]!=nullptr || node->branch[1]!=nullptr){
        uint8_t bit = get_bit(BS);
        x3f_huffnode_t *new_node = node->branch[bit];
        node = new_node;
        if (node==nullptr){
            /* TODO: Shouldn't this be treated as a fatal error? */
            x3f_printf(ERR, "Huffman coding got unexpected bit\n");
            return 0;
        }
    }
    diff = node->leaf;
    return diff;
}


static void huffman_decode_row(x3f_info_t *I, x3f_directory_entry_t *DE, int bits,
                               int row, int offset, int *minimum){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    x3f_huffman_t *HUF = ID->huffman;
    int16_t c[3] = {offset, offset, offset};
    int col;
    bit_state_t BS;

    // void pointer arithmetic
    set_bit_state(&BS, static_cast<uint8_t *>(ID->data) + HUF->row_offsets.element[row]);

    for (col=0; col<ID->columns; col++){
        int color;
        for (color=0; color<3; color++){
            uint16_t c_fix;

            c[color] += get_huffman_diff(&BS, &HUF->tree);
            if (c[color]<0){
                c_fix = 0;
                if (c[color]<*minimum)
                    *minimum = c[color];
            } else {
                c_fix = c[color];
            }

            switch (ID->type_format){
                case X3F_IMAGE_RAW_HUFFMAN_X530:
                case X3F_IMAGE_RAW_HUFFMAN_10BIT:
                    HUF->x3rgb16.data[3*(row*ID->columns + col) + color] = (uint16_t)c_fix;
                    break;
                case X3F_IMAGE_THUMB_HUFFMAN:
                    HUF->rgb8.data[3*(row*ID->columns + col) + color] = (uint8_t)c_fix;
                    break;
                default:
                    /* TODO: Shouldn't this be treated as a fatal error? */
                    x3f_printf(ERR, "Unknown huffman image type\n");
            }
        }
    }
}


static void huffman_decode(x3f_info_t *I, x3f_directory_entry_t *DE, int bits){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    int row;
    int minimum = 0;
    int offset = legacy_offset;

    x3f_printf(DEBUG, "Huffman decode with offset: %d\n", offset);
    for (row=0; row<ID->rows; row++)
        huffman_decode_row(I, DE, bits, row, offset, &minimum);

    if (auto_legacy_offset && minimum<0){
        offset = -minimum;
        x3f_printf(DEBUG, "Redo with offset: %d\n", offset);
        for (row = 0; row < ID->rows; row++)
            huffman_decode_row(I, DE, bits, row, offset, &minimum);
    }
}


static int32_t get_simple_diff(x3f_huffman_t *HUF, uint16_t index){
    if (HUF->mapping.size==0)
        return index;
    else
        return HUF->mapping.element[index];
}


static void simple_decode_row(x3f_info_t *I, x3f_directory_entry_t *DE,
                              int bits, int row, int row_stride){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    x3f_huffman_t *HUF = ID->huffman;
    // void pointer arithmetic
    uint32_t *data = (uint32_t *)(static_cast<uint8_t*>(ID->data) + row * row_stride);
    uint16_t c[3] = {0,0,0};
    int col;
    uint32_t mask = 0;

    switch (bits){
        case 8:
            mask = 0x0ff;
            break;
        case 9:
            mask = 0x1ff;
            break;
        case 10:
            mask = 0x3ff;
            break;
        case 11:
            mask = 0x7ff;
            break;
        case 12:
            mask = 0xfff;
            break;
        default:
            /* TODO: Shouldn't this be treated as a fatal error? */
            x3f_printf(ERR, "Unknown number of bits: %d\n", bits);
            mask = 0;
            break;
    }

    for (col=0; col<ID->columns; col++){
        int color;
        uint32_t val = data[col];
        for (color=0; color<3; color++){
            uint16_t c_fix;
            c[color] += get_simple_diff(HUF, (val>>(color * bits))&mask);

            switch (ID->type_format){
                case X3F_IMAGE_RAW_HUFFMAN_X530:
                case X3F_IMAGE_RAW_HUFFMAN_10BIT:
                    c_fix = (int16_t)c[color] > 0 ? c[color] : 0;
                    HUF->x3rgb16.data[3*(row*ID->columns + col) + color] = c_fix;
                    break;
                case X3F_IMAGE_THUMB_HUFFMAN:
                    c_fix = (int8_t)c[color] > 0 ? c[color] : 0;
                    HUF->rgb8.data[3*(row*ID->columns + col) + color] = c_fix;
                    break;
                default:
                    /* TODO: Shouldn't this be treated as a fatal error? */
                    x3f_printf(ERR, "Unknown huffman image type\n");
            }
        }
    }
} // simple_decode_row


static void simple_decode(x3f_info_t *I, x3f_directory_entry_t *DE,
                          int bits, int row_stride){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    int row;
    for (row=0; row<ID->rows; row++)
        simple_decode_row(I, DE, bits, row, row_stride);
}


/* --------------------------------------------------------------------- */
/* Loading the data in a directory entry                                 */
/* --------------------------------------------------------------------- */

/* First you set the offset to where to start reading the data ... */
static void read_data_set_offset(x3f_info_t *I, x3f_directory_entry_t *DE,
                                 uint32_t header_size){
    uint32_t i_off = DE->input.offset + header_size;
    fseek(I->input.file, i_off, SEEK_SET);
}


/* ... then you read the data, block for block */
static uint32_t read_data_block(void **data, x3f_info_t *I,
                                x3f_directory_entry_t *DE, uint32_t footer){
      uint32_t size = DE->input.size + DE->input.offset - ftell(I->input.file) - footer;
      *data = (void *)malloc(size);
      GETN(*data, size, I);
      return size;
}


static void x3f_load_image_verbatim(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    x3f_printf(DEBUG, "Load image verbatim\n");
    ID->data_size = read_data_block(&ID->data, I, DE, 0);
}


#if defined(_WIN32) || defined (_WIN64)
static char *utf16le_to_utf8(utf16_t *str){
    size_t osize = WideCharToMultiByte(CP_UTF8, 0, str, -1, nullptr, 0, nullptr, nullptr);
    char *buf = static_cast<char *> malloc(osize);
    WideCharToMultiByte(CP_UTF8, 0, str, -1, buf, osize, nullptr, nullptr);
    return buf;
}
#else
static char *utf16le_to_utf8(utf16_t *str){
    iconv_t ic = iconv_open("UTF-8", "UTF-16LE");
    size_t isize, osize;
    char *buf, *ibuf, *obuf;

    assert(ic != (iconv_t)-1);

    for (isize=0; str[isize]; isize++);
    isize *= 2;			/* Size in bytes */
    osize = 2*isize;		/* Worst case scenario */

    buf = static_cast<char *>(malloc(osize+1));
    ibuf = (char *)str;
    obuf = buf;

    int res = iconv(ic, &ibuf, &isize, &obuf, &osize);
    assert(res != -1);
    *obuf = 0;

    iconv_close(ic);
    return (char *) realloc(buf, obuf-buf+1);
}
#endif


static void x3f_load_property_list(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_property_list_t *PL = &DEH->data_subsection.property_list;
    int i;
    
    read_data_set_offset(I, DE, X3F_PROPERTY_LIST_HEADER_SIZE);
    GET_PROPERTY_TABLE(&PL->property_table, PL->num_properties, I);
    PL->data_size = read_data_block(&PL->data, I, DE, 0);
    
    for (i=0; i<PL->num_properties; i++){
        x3f_property_t *P = &PL->property_table.element[i];
        P->name = ((utf16_t *)PL->data + P->name_offset);
        P->value = ((utf16_t *)PL->data + P->value_offset);
        P->name_utf8 = utf16le_to_utf8(P->name);
        P->value_utf8 = utf16le_to_utf8(P->value);
    }
}


static void x3f_load_true(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    x3f_true_t *TRU = new_true(&ID->tru);
    x3f_quattro_t *Q = nullptr;
    int i;

    if (ID->type_format == X3F_IMAGE_RAW_QUATTRO ||
            ID->type_format == X3F_IMAGE_RAW_SDQ ||
            ID->type_format == X3F_IMAGE_RAW_SDQH){
        x3f_printf(DEBUG, "Load Quattro extra info\n");
        Q = new_quattro(&ID->quattro);

        for (i=0; i<TRUE_PLANES; i++){
            GET2(&Q->plane[i].columns, I);
            GET2(&Q->plane[i].rows, I);
        }

        if (Q->plane[0].rows==ID->rows/2){
            x3f_printf(DEBUG, "Quattro layout\n");
            Q->quattro_layout = 1;
        } else if (Q->plane[0].rows==ID->rows){
            x3f_printf(DEBUG, "Binned Quattro\n");
            Q->quattro_layout = 0;
        } else {
            x3f_printf(ERR, "Quattro file with unknown layer size\n");
            assert(0);
        }
    }

    x3f_printf(DEBUG, "Load TRUE\n");

    /* Read TRUE header data */
    GET2(&TRU->seed[0], I);
    GET2(&TRU->seed[1], I);
    GET2(&TRU->seed[2], I);
    GET2(&TRU->unknown, I);
    GET_TRUE_HUFF_TABLE(&TRU->table, I);

    if (ID->type_format == X3F_IMAGE_RAW_QUATTRO ||
            ID->type_format == X3F_IMAGE_RAW_SDQ ||
            ID->type_format == X3F_IMAGE_RAW_SDQH){
        x3f_printf(DEBUG, "Load Quattro extra info 2\n");
        GET4(&Q->unknown, I);
    }

    GET_TABLE_32(TRU->plane_size, TRUE_PLANES, I);

    /* Read image data */
    ID->data_size = read_data_block(&ID->data, I, DE, 0);

    /* TODO: can it be fewer than 8 bits? Maybe taken from TRU->table? */
    new_huffman_tree(&TRU->tree, 8);
    populate_true_huffman_tree(&TRU->tree, &TRU->table);

    #ifdef DBG_PRNT
    print_huffman_tree(TRU->tree.nodes, 0, 0);
    #endif

    TRU->plane_address[0] = static_cast<uint8_t *> (ID->data);
    for (i=1; i<TRUE_PLANES; i++)
        TRU->plane_address[i] = TRU->plane_address[i-1] +
                                (((TRU->plane_size.element[i-1] + 15) / 16) * 16);

    if ((ID->type_format == X3F_IMAGE_RAW_QUATTRO ||
            ID->type_format == X3F_IMAGE_RAW_SDQ ||
            ID->type_format == X3F_IMAGE_RAW_SDQH ) &&
            Q->quattro_layout){
        uint32_t columns = Q->plane[0].columns;
        uint32_t rows = Q->plane[0].rows;
        uint32_t channels = 3;
        uint32_t size = columns * rows * channels;

        TRU->x3rgb16.columns = columns;
        TRU->x3rgb16.rows = rows;
        TRU->x3rgb16.channels = channels;
        TRU->x3rgb16.row_stride = columns * channels;
        TRU->x3rgb16.buf = TRU->x3rgb16.data = static_cast<uint16_t *>(malloc(sizeof(uint16_t)*size));

        columns = Q->plane[2].columns;
        rows = Q->plane[2].rows;
        channels = 1;
        size = columns * rows * channels;

        Q->top16.columns = columns;
        Q->top16.rows = rows;
        Q->top16.channels = channels;
        Q->top16.row_stride = columns * channels;
        Q->top16.buf = Q->top16.data = static_cast<uint16_t *>(malloc(sizeof(uint16_t)*size));
    } else {
        uint32_t size = ID->columns * ID->rows * 3;
        TRU->x3rgb16.columns = ID->columns;
        TRU->x3rgb16.rows = ID->rows;
        TRU->x3rgb16.channels = 3;
        TRU->x3rgb16.row_stride = ID->columns * 3;
        TRU->x3rgb16.buf = TRU->x3rgb16.data = static_cast<uint16_t *>(malloc(sizeof(uint16_t)*size));
    }
    true_decode(I, DE);
} // x3f_load_true


static void x3f_load_huffman_compressed(x3f_info_t *I, x3f_directory_entry_t *DE,
                                        int bits, int use_map_table){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    x3f_huffman_t *HUF = ID->huffman;
    int table_size = 1<<bits;
    int row_offsets_size = ID->rows * sizeof(HUF->row_offsets.element[0]);

    x3f_printf(DEBUG, "Load huffman compressed\n");

    GET_TABLE_32(HUF->table, table_size, I);
    ID->data_size = read_data_block(&ID->data, I, DE, row_offsets_size);
    GET_TABLE_32(HUF->row_offsets, ID->rows, I);

    x3f_printf(DEBUG, "Make huffman tree ...\n");
    new_huffman_tree(&HUF->tree, bits);
    populate_huffman_tree(&HUF->tree, &HUF->table, &HUF->mapping);
    x3f_printf(DEBUG, "... DONE\n");

#ifdef DBG_PRNT
    print_huffman_tree(HUF->tree.nodes, 0, 0);
#endif
    huffman_decode(I, DE, bits);
}


static void x3f_load_huffman_not_compressed(x3f_info_t *I, x3f_directory_entry_t *DE,
                                            int bits, int use_map_table, int row_stride){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    x3f_printf(DEBUG, "Load huffman not compressed\n");
    ID->data_size = read_data_block(&ID->data, I, DE, 0);
    simple_decode(I, DE, bits, row_stride);
}


static void x3f_load_huffman(x3f_info_t *I, x3f_directory_entry_t *DE,
                             int bits, int use_map_table, int row_stride){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;
    x3f_huffman_t *HUF = new_huffman(&ID->huffman);
    uint32_t size;

    if (use_map_table){
        int table_size = 1<<bits;
        GET_TABLE_16(HUF->mapping, table_size, I);
    }

    switch (ID->type_format){
        case X3F_IMAGE_RAW_HUFFMAN_X530:
        case X3F_IMAGE_RAW_HUFFMAN_10BIT:
            size = ID->columns * ID->rows * 3;
            HUF->x3rgb16.columns = ID->columns;
            HUF->x3rgb16.rows = ID->rows;
            HUF->x3rgb16.channels = 3;
            HUF->x3rgb16.row_stride = ID->columns * 3;
            HUF->x3rgb16.buf = HUF->x3rgb16.data = static_cast<uint16_t *>(malloc(sizeof(uint16_t) * size));
            break;
        case X3F_IMAGE_THUMB_HUFFMAN:
            size = ID->columns * ID->rows * 3;
            HUF->rgb8.columns = ID->columns;
            HUF->rgb8.columns = ID->rows;
            HUF->rgb8.channels = 3;
            HUF->rgb8.row_stride = ID->columns * 3;
            HUF->rgb8.buf = HUF->rgb8.data = static_cast<uint8_t *>(malloc(sizeof(uint8_t) * size));
            break;
        default:
            /* TODO: Shouldn't this be treated as a fatal error? */
            x3f_printf(ERR, "Unknown huffman image type\n");
    }
    if (row_stride==0)
        return x3f_load_huffman_compressed(I, DE, bits, use_map_table);
    else
        return x3f_load_huffman_not_compressed(I, DE, bits, use_map_table, row_stride);
}


static void x3f_load_pixmap(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_printf(DEBUG, "Load pixmap\n");
    x3f_load_image_verbatim(I, DE);
}


static void x3f_load_jpeg(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_printf(DEBUG, "Load JPEG\n");
    x3f_load_image_verbatim(I, DE);
}


static void x3f_load_image(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_image_data_t *ID = &DEH->data_subsection.image_data;

    read_data_set_offset(I, DE, X3F_IMAGE_HEADER_SIZE);

    switch (ID->type_format){
        case X3F_IMAGE_RAW_TRUE:
        case X3F_IMAGE_RAW_MERRILL:
        case X3F_IMAGE_RAW_QUATTRO:
        case X3F_IMAGE_RAW_SDQ:
        case X3F_IMAGE_RAW_SDQH:
            x3f_load_true(I, DE);
            break;
        case X3F_IMAGE_RAW_HUFFMAN_X530:
        case X3F_IMAGE_RAW_HUFFMAN_10BIT:
            x3f_load_huffman(I, DE, 10, 1, ID->row_stride);
            break;
        case X3F_IMAGE_THUMB_PLAIN:
            x3f_load_pixmap(I, DE);
            break;
        case X3F_IMAGE_THUMB_HUFFMAN:
            x3f_load_huffman(I, DE, 8, 0, ID->row_stride);
            break;
        case X3F_IMAGE_THUMB_JPEG:
            x3f_load_jpeg(I, DE);
            break;
        default:
            /* TODO: Shouldn't this be treated as a fatal error? */
            x3f_printf(ERR, "Unknown image type\n");
    }
}


static void x3f_load_camf_decode_type2(x3f_camf_t *CAMF){
    uint32_t key = CAMF->t2.crypt_key;
    int i;

    CAMF->decoded_data_size = CAMF->data_size;
    CAMF->decoded_data = malloc(CAMF->decoded_data_size);

    for (i=0; i<CAMF->data_size; i++){
        uint8_t old, _new;
        uint32_t tmp;
        old = ((uint8_t *)CAMF->data)[i];
        key = (key * 1597 + 51749) % 244944;
        tmp = (uint32_t)(key * ((int64_t)301593171) >> 24);
        _new = (uint8_t)(old ^ (uint8_t)(((((key << 8) - tmp) >> 1) + tmp) >> 17));
        ((uint8_t *)CAMF->decoded_data)[i] = _new;
    }
}


/* NOTE: the unpacking in this code is in big respects identical to
   true_decode_one_color(). The difference is in the output you
   build. It might be possible to make some parts shared. NOTE ALSO:
   This means that the meta data is obfuscated using an image
   compression algorithm. */
static void camf_decode_type4(x3f_camf_t *CAMF){
    uint32_t seed = CAMF->t4.decode_bias;
    int row;
    uint8_t *dst;
    uint32_t dst_size = CAMF->t4.decoded_data_size;
    uint8_t *dst_end;
    bool_t odd_dst = 0;
    x3f_hufftree_t *tree = &CAMF->tree;
    bit_state_t BS;
    int32_t row_start_acc[2][2];
    uint32_t rows = CAMF->t4.block_count;
    uint32_t cols = CAMF->t4.block_size;

    CAMF->decoded_data_size = dst_size;

    CAMF->decoded_data = malloc(CAMF->decoded_data_size);
    memset(CAMF->decoded_data, 0, CAMF->decoded_data_size);

    dst = (uint8_t *)CAMF->decoded_data;
    dst_end = dst + dst_size;

    set_bit_state(&BS, CAMF->decoding_start);

    row_start_acc[0][0] = seed;
    row_start_acc[0][1] = seed;
    row_start_acc[1][0] = seed;
    row_start_acc[1][1] = seed;

    for (row=0; row<rows; row++){
        int col;
        bool_t odd_row = row&1;
        int32_t acc[2];

        /* We loop through all the columns and the rows. But the actual
           data is smaller than that, so we break the loop when reaching
           the end. */
        for (col=0; col<cols; col++){
            bool_t odd_col = col&1;
            int32_t diff = get_true_diff(&BS, tree);
            int32_t prev = col < 2 ?
            row_start_acc[odd_row][odd_col] :
            acc[odd_col];
            int32_t value = prev + diff;

            acc[odd_col] = value;
            if (col<2)
                row_start_acc[odd_row][odd_col] = value;

            switch(odd_dst){
                case 0:
                    *dst++  = (uint8_t)((value>>4)&0xff);
                    if (dst>=dst_end){
                        goto ready;
                    }
                    *dst = (uint8_t)((value<<4)&0xf0);
                    break;
                case 1:
                    *dst++ |= (uint8_t)((value>>8)&0x0f);
                    if (dst>=dst_end){
                        goto ready;
                    }
                    *dst++  = (uint8_t)((value<<0)&0xff);
                    if (dst>=dst_end){
                        goto ready;
                    }
                    break;
            }
            odd_dst = !odd_dst;
            } /* end col */
    } /* end row */
    ready:;
} // camf_decode_type4


static void x3f_load_camf_decode_type4(x3f_camf_t *CAMF){
    int i;
    uint8_t *p;
    x3f_true_huffman_element_t *element = nullptr;

    for (i=0, p=static_cast<uint8_t*>(CAMF->data); *p != 0; i++){
        /* TODO: Is this too expensive ??*/
        element = (x3f_true_huffman_element_t *)realloc(element, (i+1)*sizeof(*element));
        element[i].code_size = *p++;
        element[i].code = *p++;
    }
    CAMF->table.size = i;
    CAMF->table.element = element;

    /* TODO: where does the values 28 and 32 come from? */
    constexpr int CAMF_T4_DATA_SIZE_OFFSET = 28;
    constexpr int CAMF_T4_DATA_OFFSET = 32;
    // void pointer arithmetic
    CAMF->decoding_size = *(uint32_t *)(static_cast<uint8_t*>(CAMF->data) + CAMF_T4_DATA_SIZE_OFFSET);
    CAMF->decoding_start = (uint8_t *)CAMF->data + CAMF_T4_DATA_OFFSET;

    /* TODO: can it be fewer than 8 bits? Maybe taken from TRU->table? */
    new_huffman_tree(&CAMF->tree, 8);
    populate_true_huffman_tree(&CAMF->tree, &CAMF->table);

#ifdef DBG_PRNT
    print_huffman_tree(CAMF->tree.nodes, 0, 0);
#endif

    camf_decode_type4(CAMF);
} // x3f_load_camf_decode_type4


static void camf_decode_type5(x3f_camf_t *CAMF){
    int32_t acc = CAMF->t5.decode_bias;
    uint8_t *dst;
    x3f_hufftree_t *tree = &CAMF->tree;
    bit_state_t BS;
    int32_t i;

    CAMF->decoded_data_size = CAMF->t5.decoded_data_size;
    CAMF->decoded_data = (void *) malloc(CAMF->decoded_data_size);
    dst = (uint8_t *)CAMF->decoded_data;
    set_bit_state(&BS, CAMF->decoding_start);

    for (i=0; i<CAMF->decoded_data_size; i++){
        int32_t diff = get_true_diff(&BS, tree);
        acc = acc + diff;
        *dst++ = (uint8_t)(acc & 0xff);
    }
} // x3f_load_camf_decode_type4


static void x3f_load_camf_decode_type5(x3f_camf_t *CAMF){
    int i;
    uint8_t *p;
    x3f_true_huffman_element_t *element = nullptr;

    for (i=0, p=static_cast<uint8_t*>(CAMF->data); *p != 0; i++){
        /* TODO: Is this too expensive ??*/
        element = (x3f_true_huffman_element_t *)realloc(element, (i+1)*sizeof(*element));
        element[i].code_size = *p++;
        element[i].code = *p++;
    }
    CAMF->table.size = i;
    CAMF->table.element = element;

    /* TODO: where does the values 28 and 32 come from? */
    constexpr int CAMF_T5_DATA_SIZE_OFFSET = 28;
    constexpr int CAMF_T5_DATA_OFFSET = 32;
    // void pointer arithmetic
    CAMF->decoding_size = *(uint32_t *)(static_cast<uint8_t*>(CAMF->data) + CAMF_T5_DATA_SIZE_OFFSET);
    CAMF->decoding_start = (uint8_t *)CAMF->data + CAMF_T5_DATA_OFFSET;

    /* TODO: can it be fewer than 8 bits? Maybe taken from TRU->table? */
    new_huffman_tree(&CAMF->tree, 8);
    populate_true_huffman_tree(&CAMF->tree, &CAMF->table);

#ifdef DBG_PRNT
    print_huffman_tree(CAMF->tree.nodes, 0, 0);
#endif

    camf_decode_type5(CAMF);
} // x3f_load_camf_decode_type5


static void x3f_setup_camf_text_entry(camf_entry_t *entry){
    entry->text_size = *(uint32_t *)entry->value_address;
    // void pointer arithmetic
    entry->text = static_cast<char *>(entry->value_address) + 4;
}


static void x3f_setup_camf_property_entry(camf_entry_t *entry){
    int i;
    uint8_t *e = static_cast<uint8_t*>(entry->entry);
    uint8_t *v = static_cast<uint8_t*>(entry->value_address);
    uint32_t num = entry->property_num = *(uint32_t *)v;
    uint32_t off = *(uint32_t *)(v + 4);

    entry->property_name = (char **)malloc(num*sizeof(uint8_t*));
    entry->property_value = (uint8_t **)malloc(num*sizeof(uint8_t*));

    for (i=0; i<num; i++){
        uint32_t name_off = off + *(uint32_t *)(v + 8 + 8*i);
        uint32_t value_off = off + *(uint32_t *)(v + 8 + 8*i + 4);
        entry->property_name[i] = (char *)(e + name_off);
        entry->property_value[i] = e + value_off;
    }
}


static void set_matrix_element_info(uint32_t _type, uint32_t *size, matrix_type_t *decoded_type){
    switch (_type){
        case 0:
            *size = 2;
            *decoded_type = M_INT; /* known to be true */
            break;
        case 1:
            *size = 4;
            *decoded_type = M_UINT; /* TODO: unknown ???? */
            break;
        case 2:
            *size = 4;
            *decoded_type = M_UINT; /* TODO: unknown ???? */
            break;
        case 3:
            *size = 4;
            *decoded_type = M_FLOAT; /* known to be true */
            break;
        case 5:
            *size = 1;
            *decoded_type = M_UINT; /* TODO: unknown ???? */
            break;
        case 6:
            *size = 2;
            *decoded_type = M_UINT; /* TODO: unknown ???? */
            break;
        default:
            x3f_printf(ERR, "Unknown matrix type (%ud)\n", _type);
            assert(0);
    }
}


static void get_matrix_copy(camf_entry_t *entry){
    uint32_t element_size = entry->matrix_element_size;
    uint32_t elements = entry->matrix_elements;
    int i, size = (entry->matrix_decoded_type==M_FLOAT ?
                  sizeof(double) : sizeof(uint32_t)) * elements;
    entry->matrix_decoded = (void *) malloc(size);

    switch (element_size){
        case 4:
            switch (entry->matrix_decoded_type){
                case M_INT:
                case M_UINT:
                    memcpy(entry->matrix_decoded, entry->matrix_data, size);
                    break;
                case M_FLOAT:
                    for (i=0; i<elements; i++)
                        ((double *)entry->matrix_decoded)[i] =
                        (double)((float *)entry->matrix_data)[i];
                    break;
                default:
                    x3f_printf(ERR, "Invalid matrix element type of size 4\n");
                    assert(0);
            }
            break;
        case 2:
            switch (entry->matrix_decoded_type){
                case M_INT:
                    for (i=0; i<elements; i++)
                        ((int32_t *)entry->matrix_decoded)[i] =
                        (int32_t)((int16_t *)entry->matrix_data)[i];
                    break;
                case M_UINT:
                    for (i=0; i<elements; i++)
                        ((uint32_t *)entry->matrix_decoded)[i] =
                        (uint32_t)((uint16_t *)entry->matrix_data)[i];
                    break;
                default:
                    x3f_printf(ERR, "Invalid matrix element type of size 2\n");
                    assert(0);
            }
            break;
        case 1:
            switch (entry->matrix_decoded_type){
                case M_INT:
                    for (i=0; i<elements; i++)
                        ((int32_t *)entry->matrix_decoded)[i] =
                        (int32_t)((int8_t *)entry->matrix_data)[i];
                    break;
                case M_UINT:
                    for (i=0; i<elements; i++)
                        ((uint32_t *)entry->matrix_decoded)[i] =
                        (uint32_t)((uint8_t *)entry->matrix_data)[i];
                    break;
                default:
                    x3f_printf(ERR, "Invalid matrix element type of size 1\n");
                    assert(0);
            }
            break;
        default:
            x3f_printf(ERR, "Unknown size %d\n", element_size);
            assert(0);
    }
} // get_matrix_copy


static void x3f_setup_camf_matrix_entry(camf_entry_t *entry){
    int i;
    int totalsize = 1;

    uint8_t *e = static_cast<uint8_t*>(entry->entry);
    uint8_t *v = static_cast<uint8_t*>(entry->value_address);
    uint32_t _type = entry->matrix_type = *(uint32_t *)(v + 0);
    uint32_t dim = entry->matrix_dim = *(uint32_t *)(v + 4);
    uint32_t off = entry->matrix_data_off = *(uint32_t *)(v + 8);
    camf_dim_entry_t *dentry = entry->matrix_dim_entry =
                    (camf_dim_entry_t*)malloc(dim*sizeof(camf_dim_entry_t));

    for (i=0; i<dim; i++){
        uint32_t size = dentry[i].size = *(uint32_t *)(v + 12 + 12*i + 0);
        dentry[i].name_offset = *(uint32_t *)(v + 12 + 12*i + 4);
        dentry[i].n = *(uint32_t *)(v + 12 + 12*i + 8);
        dentry[i].name = (char *)(e + dentry[i].name_offset);

        if (dentry[i].n != i){
            /* TODO: is something needed to be made in this case */
            x3f_printf(DEBUG, "Matrix entry for %s/%s is out of order "
                    "(index/%d != order/%d)\n",
                    entry->name_address, dentry[i].name, dentry[i].n, i);
        }
        totalsize *= size;
    }

    set_matrix_element_info(_type, &entry->matrix_element_size,
                            &entry->matrix_decoded_type);
    entry->matrix_data = (void *)(e + off);
    entry->matrix_elements = totalsize;
    entry->matrix_used_space = entry->entry_size - off;

    /* This estimate only works for matrices above a certain size */
    entry->matrix_estimated_element_size = entry->matrix_used_space / totalsize;
    get_matrix_copy(entry);
} // x3f_setup_camf_matrix_entry


static void x3f_setup_camf_entries(x3f_camf_t *CAMF){
    uint8_t *p = (uint8_t *)CAMF->decoded_data;
    uint8_t *end = p + CAMF->decoded_data_size;
    camf_entry_t *entry = nullptr;
    int i;

    x3f_printf(DEBUG, "SETUP CAMF ENTRIES\n");

    for (i=0; p < end; i++){
        uint32_t *p4 = (uint32_t *)p;

        switch (*p4){
            case X3F_CMbP:
            case X3F_CMbT:
            case X3F_CMbM:
                break;
            default:
                x3f_printf(ERR, "Unknown CAMF entry %x @ %p\n", *p4, p4);
                x3f_printf(ERR, "  start = %p end = %p\n", CAMF->decoded_data, end);
                x3f_printf(ERR, "  left = %ld\n", (long)(end - p));
                x3f_printf(ERR, "Stop parsing CAMF\n");
                /* TODO: Shouldn't this be treated as a fatal error? */
                goto stop;
        }

        /* TODO: lots of realloc - may be inefficient */
        entry = (camf_entry_t *)realloc(entry, (i+1) * sizeof(camf_entry_t));

        /* Pointer */
        entry[i].entry = p;

        /* Header */
        entry[i].id = *p4++;
        entry[i].version = *p4++;
        entry[i].entry_size = *p4++;
        entry[i].name_offset = *p4++;
        entry[i].value_offset = *p4++;

        /* Compute adresses and sizes */
        entry[i].name_address = (char *)(p + entry[i].name_offset);
        entry[i].value_address = p + entry[i].value_offset;
        entry[i].name_size = entry[i].value_offset - entry[i].name_offset;
        entry[i].value_size = entry[i].entry_size - entry[i].value_offset;

        entry[i].text_size = 0;
        entry[i].text = nullptr;
        entry[i].property_num = 0;
        entry[i].property_name = nullptr;
        entry[i].property_value = nullptr;
        entry[i].matrix_type = 0;
        entry[i].matrix_dim = 0;
        entry[i].matrix_data_off = 0;
        entry[i].matrix_data = nullptr;
        entry[i].matrix_dim_entry = nullptr;

        entry[i].matrix_decoded = nullptr;

        switch (entry[i].id){
            case X3F_CMbP:
                x3f_setup_camf_property_entry(&entry[i]);
                break;
            case X3F_CMbT:
                x3f_setup_camf_text_entry(&entry[i]);
                break;
            case X3F_CMbM:
                x3f_setup_camf_matrix_entry(&entry[i]);
                break;
        }
        p += entry[i].entry_size;
    }

    stop:
    CAMF->entry_table.size = i;
    CAMF->entry_table.element = entry;
    x3f_printf(DEBUG, "SETUP CAMF ENTRIES (READY) Found %d entries\n", i);
} // x3f_setup_camf_entries


static void x3f_load_camf(x3f_info_t *I, x3f_directory_entry_t *DE){
    x3f_directory_entry_header_t *DEH = &DE->header;
    x3f_camf_t *CAMF = &DEH->data_subsection.camf;

    x3f_printf(DEBUG, "Loading CAMF of type %d\n", CAMF->_type);
    read_data_set_offset(I, DE, X3F_CAMF_HEADER_SIZE);
    CAMF->data_size = read_data_block(&CAMF->data, I, DE, 0);

    switch (CAMF->_type){
        case 2:			/* Older SD9-SD14 */
            x3f_load_camf_decode_type2(CAMF);
            break;
        case 4:			/* TRUE ... Merrill */
            x3f_load_camf_decode_type4(CAMF);
            break;
        case 5:			/* Quattro ... */
            x3f_load_camf_decode_type5(CAMF);
            break;
        default:
            /* TODO: Shouldn't this be treated as a fatal error? */
            x3f_printf(ERR, "Unknown CAMF type\n");
    }

    if (CAMF->decoded_data != nullptr)
        x3f_setup_camf_entries(CAMF);
    else
        /* TODO: Shouldn't this be treated as a fatal error? */
        x3f_printf(ERR, "No decoded CAMF data\n");
} // x3f_load_camf


/* extern */ 
x3f_return_t x3f_load_data(x3f_t *x3f, x3f_directory_entry_t *DE){
    x3f_info_t *I = &x3f->info;

    if (DE==nullptr)
        return X3F_ARGUMENT_ERROR;

    switch (DE->header.identifier){
        case X3F_SECp:
            x3f_load_property_list(I, DE);
            break;
        case X3F_SECi:
            x3f_load_image(I, DE);
            break;
        case X3F_SECc:
            x3f_load_camf(I, DE);
            break;
        default:
            x3f_printf(ERR, "Unknown directory entry type\n");
            return X3F_INTERNAL_ERROR;
    }
    return X3F_OK;
}


/* extern */ 
x3f_return_t x3f_load_image_block(x3f_t *x3f, x3f_directory_entry_t *DE){
    x3f_info_t *I = &x3f->info;

    if (DE==nullptr)
        return X3F_ARGUMENT_ERROR;
    x3f_printf(DEBUG, "Load image block\n");

    switch (DE->header.identifier){
        case X3F_SECi:
            read_data_set_offset(I, DE, X3F_IMAGE_HEADER_SIZE);
            x3f_load_image_verbatim(I, DE);
            break;
        default:
            x3f_printf(ERR, "Unknown image directory entry type\n");
            return X3F_INTERNAL_ERROR;
    }
    return X3F_OK;
}


/* extern */ 
const char *x3f_err(x3f_return_t err){
    switch (err){
        case X3F_OK:             return "ok";
        case X3F_ARGUMENT_ERROR: return "argument error";
        case X3F_INFILE_ERROR:   return "infile error";
        case X3F_OUTFILE_ERROR:  return "outfile error";
        case X3F_INTERNAL_ERROR: return "internal error";
        default:                 return "undefined error";
    }
}

/* --------------------------------------------------------------------- */
/* The End                                                               */
/* --------------------------------------------------------------------- */


} // namespace x3ftools
