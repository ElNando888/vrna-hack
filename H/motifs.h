#ifndef __VIENNA_RNA_PACKAGE_MOTIFS_H__
#define __VIENNA_RNA_PACKAGE_MOTIFS_H__

#include "data_structures.h"

typedef struct _ligand {
  /* static fields */
  const char *name;
  FLT_OR_DBL Kd;
  /* runtime */
  FLT_OR_DBL conc;        /* in micromolars */ 
  int        deltaG;      /* in dcal/mol    */
} ligand;

typedef struct _motif {
  /* static fields */
  const char *name;
  unsigned   num_segments;
  const char **segment;   /* for matching against the sequence */
  int        *s_ofs;      /* offset from th matched substring */
  int        *s_len;      /* for comparison in the callbacks */
  int        intrinsic;   /* in dcal/mol, estimated deviation from the model */
  int        lig_index;   /* which ligand, if any */
  /* runtime */
  int        **occur;
  int        deltaG;      /* in dcal/mol */
} motif;


ligand* get_ligands(void);
int set_ligand(ligand* lig_db, const char* lname, FLT_OR_DBL concentration);

motif* get_motifs(void);
void reset_motifs(motif* mdb);
void detect_motifs(const char *sequence, motif* mdb, ligand* lig_db);

void setup_motifs_for_params(paramT* P, motif* mdb);
void setup_motifs_for_pf_params(pf_paramT* P, motif* mdb);

#endif
