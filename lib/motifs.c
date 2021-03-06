/*
                  Motifs, aptamers and ligands
                  
                  F Portela
                  Vienna RNA "hack"
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>  
#include <ctype.h> 
#include <unistd.h>
#include <string.h>
#include "fold_vars.h"
#include "motifs.h"

#define PUBLIC
#define PRIVATE static

typedef enum {false, true} bool;

/*--------------------------------------------------------------------------*/

PRIVATE ligand known_ligands[] = {
  {"FMN", 0.},
  {"TEP", 0.},
  {"MS2", 0.},
  {"ARG", 0.},
  {"TRP", 0.},
  {"ATP", 0.}
};

int num_ligands = sizeof(known_ligands) / sizeof(ligand);

enum {
  _FMN,
  _TEP,
  _MS2,
  _ARG,
  _TRP,
  _ATP,
  _NONE = -1
};

/* TODO: retrieve these from a file */
PRIVATE motif  known_motifs[] = {
  // .intrinsic estimated from R91 experimental results
  {"FMN aptamer", 2, (const char *[]){"AGGAUAU","AGAAGG"}, (int[]){0,1}, (int[]){6,5}, -60, 0, _FMN, 0.3, NULL, 0},
  {"FMN aptamer", 2, (const char *[]){"AGGAUAN","NGAAGG"}, (int[]){0,1}, (int[]){6,5}, -60, 0, _FMN, 3.0, NULL, 0},
  // .intrinsic estimated from R92 experimental results
  {"TEP aptamer", 2, (const char *[]){"CCUUGGCAG","AUACCA"}, (int[]){0,0}, (int[]){9,6}, 0, -440, _TEP, 0.32, NULL, 0},
  // MS2
  {"MS2 hairpin", 2, (const char *[]){"ACAUGAGGAUCACCCAUGU",NULL}, (int[]){5,14}, (int[]){1,0}, 0, 0, _MS2, 0.00256, NULL, 0},
  // FIXME: wrong Kd
  {"ARG aptamer", 2, (const char *[]){"GAAGGAGCG","CAGGUAGGUCGC"}, (int[]){1,1}, (int[]){7,10}, 0, -160, _ARG, 0.32, NULL, 0},
  // FIXME: wrong Kd
  {"TRP hairpin", 2, (const char *[]){"CGAGGACCGGUACGGCCGCCACUCG",NULL}, (int[]){4,18}, (int[]){3,3}, 0, 0, _TRP, 0.32, NULL, 0},
  {"TRP aptamer", 2, (const char *[]){"CGGCCGCCACU","AGGACCGGG"}, (int[]){1,2}, (int[]){8,6}, 0, -400, _TRP, 0.32, NULL, 0},
  // FIXME: wrong Kd
  {"ATP aptamer", 2, (const char *[]){"GGGAAGAAACUGC","GGC"}, (int[]){1,1}, (int[]){11,1}, 0, 0, _ATP, 0.32, NULL, 0},
  // meant for tests, no experimental data (so far) to support the .intrinsic value below
  // {"Sarcin-ricin (example)", 2, (const char*[]){"CCAGUA","GAACA"}, (int[]){0,0}, (int[]){6,5}, -250, _NONE, 0., NULL, 0}
};

int num_motifs = sizeof(known_motifs) / sizeof(motif);

/*--------------------------------------------------------------------------*/

#if 1
#define _RT (-(temperature+K0)*GASCONST/1000.)
#else
#define _RT -0.6
#endif

PUBLIC int set_ligand(ligand* lig_db, const char* lname, FLT_OR_DBL concentration)
{
  int i;
  for (i=0; i < num_ligands; i++) {
    ligand* kli = lig_db+i;
    if (strcmp(lname, kli->name)!=0) continue;
    kli->conc = concentration;
    return i;
  }
  return -1;
}


PUBLIC ligand* get_ligands(void)
{
  ligand* p = (ligand*) calloc(num_ligands, sizeof(ligand));
  if (p) memmove(p, known_ligands, num_ligands * sizeof(ligand));
  return p;
}


PUBLIC motif* get_motifs(void)
{
  motif* p = (motif*) calloc(num_motifs, sizeof(motif));
  if (p) memmove(p, known_motifs, num_motifs * sizeof(motif));
  return p;
}


PUBLIC void reset_motifs(motif* mdb)
{
  int i,k;
  for (i = 0; i < num_motifs; i++) {
    motif* kmi = mdb+i;
    if (kmi->occur == NULL) continue;
    for (k = 0; k < kmi->num_segments; k++) {
      if (kmi->occur[k]) free(kmi->occur[k]);
    }
    free(kmi->occur);
    kmi->occur = NULL;
    kmi->deltaG = 0;
  }
}


PRIVATE int* add_list(int* list, int val)
{
  int size = 1+list[0];
  list = realloc(list, (size+1)*sizeof(int));
  list[0]++;
  list[list[0]] = val;
  return list;
}


PRIVATE int in_list(int* list, int val)
{
  int i;
  for(i=1; i<=list[0]; i++) if (list[i]==val) return i;
  return 0;
}


PRIVATE bool iupac_cmp(const char nt, const char mask)
{
  char* p = NULL;
  switch(nt) {
  case 'A': p = strchr("ARMWVHDN", mask); break;
  case 'U': p = strchr("UYKWBHDN", mask); break;
  case 'G': p = strchr("GRKSBVDN", mask); break;
  case 'C': p = strchr("CYMSBVHN", mask); break;
  }
  return (p != NULL);
}

/* we may want to add this one to utils.c some day... 
 */
PRIVATE const char* iupac_match(const char* seq, const char* iupac_mask)
{
  const char *a, *b;

  b = iupac_mask;
  if (*b == 0) return seq;
  for ( ; *seq != 0; seq++) {
    if (!iupac_cmp(*seq, *b)) continue;
    a = seq;
    while (1) {
      if (*b == 0) return seq;
      if (!iupac_cmp(*a++, *b++)) break;
    }
    b = iupac_mask;
  }
  return NULL;
}


PUBLIC void detect_motifs(const char *sequence, motif* mdb, ligand* lig_db)
{
  int i,j,k;
  
  reset_motifs(mdb);
  
  for (i = 0; i < num_motifs; i++) {
    motif* kmi = mdb+i;
    for (j = 0; j < kmi->num_segments; j++) {
      const char* p;
      const char* needle = kmi->segment[j];
      if (needle == NULL) continue;
      for (p = iupac_match(sequence, needle); p; p = iupac_match(p+1, needle)) {
        int ofs = p - sequence;
        if (!ofs) continue;
        if (kmi->occur == NULL) {
          kmi->occur = (int**)calloc(kmi->num_segments, sizeof(int*));
          for (k = 0; k < kmi->num_segments; k++) {
            kmi->occur[k] = (int*)calloc(1, sizeof(int));
          }
        }
        if (0) fprintf(stderr,"%s[%d] found at %d\n", kmi->name, j, ofs+1);
        kmi->occur[j] = add_list(kmi->occur[j], ofs+1 + kmi->s_ofs[j]);
        if (kmi->lig_index >= 0 && lig_db[kmi->lig_index].conc > kmi->Kd) {
          kmi->deltaG = 100. * (_RT * log( lig_db[kmi->lig_index].conc / kmi->Kd));
        }

        // special case
        if (j == 0 && kmi->num_segments == 2 && kmi->segment[1] == NULL) {
          kmi->occur[1] = add_list(kmi->occur[1], ofs+1 + kmi->s_ofs[1]);
        }
      }
    }
  }
}


PRIVATE void std_eilcb(int* fe, int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, int ii, int qq, paramT *P)
{
  int i;
  motif* mdb;
  if (P && P->userdata) {
    mdb = (motif*)P->userdata;
  } else {
    return;
  }
  for (i = 0; i < num_motifs; i++) {
    motif* kmi = mdb+i;
    if (kmi->num_segments != 2) continue;
    if (kmi->occur == NULL) continue;

    if (kmi->s_len[0]==n1 && kmi->s_len[1]==n2) {
      if (0) fprintf(stderr,"n1-n2 %d,%d\n", ii, qq);
      if (in_list(kmi->occur[0], ii) && in_list(kmi->occur[1], qq)) {
        (*fe) += kmi->intrinsic;
        if (kmi->lig_index >= 0) {
          (*fe) += kmi->deltaG;
          (*fe) += kmi->correction;
        }
        if (0) fprintf(stderr, "%s at %d+%d\n", kmi->name, ii, qq);
        return;
      }
    }
    if (kmi->s_len[0]==n2 && kmi->s_len[1]==n1) {
      if (0) fprintf(stderr,"n2-n1 %d,%d\n", qq, ii);
      if (in_list(kmi->occur[0], qq) && in_list(kmi->occur[1], ii)) {
        (*fe) += kmi->intrinsic;
        if (kmi->lig_index >= 0) {
          (*fe) += kmi->deltaG;
          (*fe) += kmi->correction;
        }
        if (0) fprintf(stderr, "%s at %d+%d\n", kmi->name, qq, ii);
        return;
      }
    }
  }
}


PRIVATE void std_eeilcb(double* fe, int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, int ii, int qq, pf_paramT *P)
{
  int i;
  motif* mdb;
  if (P && P->userdata) {
    mdb = (motif*)P->userdata;
  } else {
    return;
  }
  for (i = 0; i < num_motifs; i++) {
    motif* kmi = mdb+i;
    if (kmi->num_segments != 2) continue;
    if (kmi->occur == NULL) continue;

    if (kmi->s_len[0]==u1 && kmi->s_len[1]==u2) {
      if (in_list(kmi->occur[0], ii) && in_list(kmi->occur[1], qq)) {
        (*fe) *= exp(-kmi->intrinsic*10./(P->kT));
        if (kmi->lig_index >= 0) {
          (*fe) *= exp(-kmi->deltaG*10./(P->kT));
          (*fe) *= exp(-kmi->correction*10./(P->kT));
        }
        return;
      }
    }
    if (kmi->s_len[0]==u2 && kmi->s_len[1]==u1) {
      if (in_list(kmi->occur[0], qq) && in_list(kmi->occur[1], ii)) {
        (*fe) *= exp(-kmi->intrinsic*10./(P->kT));
        if (kmi->lig_index >= 0) {
          (*fe) *= exp(-kmi->deltaG*10./(P->kT));
          (*fe) *= exp(-kmi->correction*10./(P->kT));
        }
        return;
      }
    }
  }
}


PUBLIC void setup_motifs_for_params(paramT* P, motif* mdb)
{
  if (P) {
    P->userdata = mdb;
    P->eilcb = std_eilcb;
  }
}


PUBLIC void setup_motifs_for_pf_params(pf_paramT* P, motif* mdb)
{
  if (P) {
    P->userdata = mdb;
    P->eeilcb = std_eeilcb;
  }
}

