#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Datatypes
typedef enum { false, true } bool;

// The atomic properties
typedef struct {
  int index;
  char* species; 
  char* label;
  double position[3];
} MagresAtom;

void magres_atom_dealloc(MagresAtom *atom) {
  if(atom->species != NULL) {
    free(atom->species);
    atom->species = NULL;
  }
  
  if(atom->label != NULL) {
    free(atom->label);
    atom->label = NULL;
  }
}

typedef struct {
  double lattice[3][3];
} MagresLattice;

typedef struct {
  char* symmetry_string;
} MagresSymmetry;

void magres_symmetry_dealloc(MagresSymmetry *symmetry) {
  if(symmetry->symmetry_string != NULL) {
    free(symmetry->symmetry_string);
    symmetry->symmetry_string = NULL;
  }
}

// The magnetic resonance properties
typedef struct {
  MagresAtom *atom;
  double sigma[3][3];
} MagresMs;

typedef struct {
  MagresAtom* atom1;
  MagresAtom* atom2;
  double K[3][3];
} MagresIsc;

typedef struct {
  MagresAtom* atom;
  double V[3][3];
} MagresEfg;

// The whole file
typedef struct {
  MagresLattice *lattice;

  int num_atoms;
  MagresAtom **atoms;

  int num_symmetries;
  MagresSymmetry **symmetries;

  int num_ms;
  MagresMs** ms;

  int num_isc;
  MagresIsc** isc;

  int num_efg;
  MagresEfg** efg;
} MagresFile;

void magres_file_init(MagresFile *magres_file) {
  magres_file->lattice = NULL;
  
  magres_file->num_atoms = 0;
  magres_file->atoms = NULL;

  magres_file->num_symmetries = 0;
  magres_file->symmetries = NULL;

  magres_file->num_ms = 0;
  magres_file->ms = NULL;

  magres_file->num_isc = 0;
  magres_file->isc = NULL;

  magres_file->num_efg = 0;
  magres_file->efg = NULL;
}

void magres_file_dealloc(MagresFile *magres_file) {
  int i;

  if(magres_file->atoms != NULL) {
    for(i = 0; i<magres_file->num_atoms; ++i) {
      magres_atom_dealloc(magres_file->atoms[i]);
      free(magres_file->atoms[i]);
    }

    free(magres_file->atoms);
  }
  
  if(magres_file->symmetries != NULL) {
    for(i = 0; i<magres_file->num_symmetries; ++i) {
      magres_symmetry_dealloc(magres_file->symmetries[i]);
      free(magres_file->symmetries[i]);
    }

    free(magres_file->symmetries);
  }

  if(magres_file->lattice != NULL) {
    free(magres_file->lattice);
    magres_file->lattice = NULL;
  }

  if(magres_file->isc != NULL) {
    for(i = 0; i<magres_file->num_isc; ++i) {
      if(magres_file->isc[i] != NULL) {
        free(magres_file->isc[i]);
        magres_file->isc[i] = 0;
      }
    }

    free(magres_file->isc);
  }

  if(magres_file->efg != NULL) {
    for(i = 0; i<magres_file->num_efg; ++i) {
      if(magres_file->efg[i] != NULL) {
        free(magres_file->efg[i]);
        magres_file->efg[i] = 0;
      }
    }

    free(magres_file->efg);
  }

  if(magres_file->ms != NULL) {
    for(i = 0; i<magres_file->num_ms; ++i) {
      if(magres_file->ms[i] != NULL) {
        free(magres_file->ms[i]);
        magres_file->ms[i] = 0;
      }
    }

    free(magres_file->ms);
  }
}

// Generic line tokenizer
typedef struct {
  char **cols;
  int num_cols;
} MagresLine;

void magres_line_dealloc(MagresLine* line) {
  int i = 0;
  for(; i<line->num_cols; ++i) {
    free(line->cols[i]);
    line->cols[i] = NULL;
  }

  free(line->cols);
  line->cols = NULL;
}

char* str_trim(char *str) {
  while(*str == ' ') {
    str++;
  }

  char *p = &str[strlen(str) - 1];

  while(*p == ' ') {
    *p = 0;
    p--;
  }

  return str;
}

MagresLine *magres_tokenize(char *str_line) {
  MagresLine *line = malloc(sizeof(MagresLine));

  char *str_line_trimmed = str_trim(str_line);
  
  // First, work out how many columns there are
  line->num_cols = 0;

  bool whitespace = false;
  int line_length = strlen(str_line_trimmed);

  int i = 0;
  for(; i<line_length; ++i) {
    if(str_line_trimmed[i] == ' ' && !whitespace) {
      whitespace = true;
      line->num_cols++;
    } else if(str_line_trimmed[i] != ' ' && whitespace) {
      whitespace = false;
    }
  }

  if(strlen(str_line_trimmed) != 0) {
    line->num_cols += 1;
  }
  
  // Now copy all the tokens into the cols array
  line->cols = malloc(sizeof(char*) * line->num_cols);

  i = 0;
  int j = 0;
  int token_length = 0;
  int current_col = 0;

  for(; i<=line_length; ++i) {
    if((str_line_trimmed[i] == ' ' && !whitespace) || i == line_length) { // start of whitespace
      token_length = i - j;

      line->cols[current_col] = malloc(sizeof(char) * (token_length + 1));
      memcpy(line->cols[current_col], &str_line_trimmed[j], token_length);
      line->cols[current_col][token_length] = 0;
     
      if(token_length != 0) {
        current_col += 1;
      }

      whitespace = true;
    } else if(str_line_trimmed[i] != ' ' && whitespace) { // end of whitespace
      j = i;
      whitespace = false;
    }
  }

  return line;
}

// [atoms] block parser
MagresAtom* magres_parse_atom(MagresLine *line) {
  if(line->num_cols != 7) {
    fprintf(stderr, "atom: Wrong number of columns, %d\n", line->num_cols);
    return NULL;
  }

  MagresAtom *atom = malloc(sizeof(MagresAtom));
  
  int species_len = strlen(line->cols[1]);

  atom->species = malloc(sizeof(char) * (species_len + 1));
  memcpy(atom->species, line->cols[1], species_len);
  atom->species[species_len] = 0;

  int label_len = strlen(line->cols[2]);

  atom->label = malloc(sizeof(char) * (label_len + 1));
  memcpy(atom->label, line->cols[2], label_len);
  atom->label[label_len] = 0;

  atom->index = atoi(line->cols[3]);

  sscanf(line->cols[4], "%lg", &atom->position[0]);
  sscanf(line->cols[5], "%lg", &atom->position[1]);
  sscanf(line->cols[6], "%lg", &atom->position[2]);

  return atom;
}

MagresSymmetry* magres_parse_symmetry(MagresLine *line) {
  if(line->num_cols != 2) {
    fprintf(stderr, "symmetry: Wrong number of columns, %d\n", line->num_cols);
    return NULL;
  }

  MagresSymmetry *symmetry = malloc(sizeof(MagresSymmetry));
  int str_len = strlen(line->cols[1]);

  symmetry->symmetry_string = malloc(sizeof(char) * (str_len + 1));
  memcpy(symmetry->symmetry_string, line->cols[1], str_len);
  symmetry->symmetry_string[str_len] = 0;

  return symmetry;
}

MagresLattice *magres_parse_lattice(MagresLine *line) {
  if(line->num_cols != 10) {
    fprintf(stderr, "lattice: Wrong number of columns, %d\n", line->num_cols);
    return NULL;
  }

  MagresLattice *lattice = malloc(sizeof(MagresLattice));
  
  sscanf(line->cols[1], "%lg", &lattice->lattice[0][0]);
  sscanf(line->cols[2], "%lg", &lattice->lattice[0][1]);
  sscanf(line->cols[3], "%lg", &lattice->lattice[0][2]);
  
  sscanf(line->cols[4], "%lg", &lattice->lattice[1][0]);
  sscanf(line->cols[5], "%lg", &lattice->lattice[1][1]);
  sscanf(line->cols[6], "%lg", &lattice->lattice[1][2]);
  
  sscanf(line->cols[7], "%lg", &lattice->lattice[2][0]);
  sscanf(line->cols[8], "%lg", &lattice->lattice[2][1]);
  sscanf(line->cols[9], "%lg", &lattice->lattice[2][2]);

  return lattice;
}

bool magres_parse_atoms(MagresLine **lines, int num_lines, MagresFile *magres_file) {
  // Count how many of each line type we have and allocate space in the MagresFile
  int num_atoms = 0, num_symmetries = 0, num_lattice = 0;
  int i_line = 0;
  MagresLine *line = NULL;

  for(i_line = 0; i_line < num_lines; ++i_line) {
    line = lines[i_line];

    if(line->num_cols != 0) {
      if(strcmp(line->cols[0], "atom") == 0) {
        num_atoms += 1;
      } else if(strcmp(line->cols[0], "symmetry") == 0) {
        num_symmetries += 1;
      } else if(strcmp(line->cols[0], "lattice") == 0) {
        num_lattice += 1;
      }
    }
  }

  if(num_atoms > 0) {
    magres_file->atoms = malloc(sizeof(MagresAtom*) * num_atoms);
    magres_file->num_atoms = num_atoms;
  } else {
    magres_file->atoms = NULL;
    magres_file->num_atoms = 0;
  }

  if(num_symmetries > 0) {
    magres_file->symmetries = malloc(sizeof(MagresSymmetry*) * num_symmetries);
    magres_file->num_symmetries = num_symmetries;
  } else {
    magres_file->symmetries = NULL;
    magres_file->num_symmetries = 0;
  }

  if(num_lattice > 0) {
    magres_file->lattice = malloc(sizeof(MagresLattice) * num_symmetries);
  } else {
    magres_file->lattice = NULL;
  }

  // Parse the lines into different types
  int i_atom=0, i_symmetry=0;

  MagresAtom *current_atom = NULL;
  MagresSymmetry *current_symmetry = NULL;

  for(i_line = 0; i_line < num_lines; ++i_line) {
    line = lines[i_line];

    if(line->num_cols != 0) {
      if(strcmp(line->cols[0], "atom") == 0) {
        current_atom = magres_parse_atom(line);
        
        if(current_atom == NULL) {
          return false;
        }

        magres_file->atoms[i_atom++] = current_atom;
      } else if(strcmp(line->cols[0], "lattice") == 0) {
        magres_file->lattice = magres_parse_lattice(line);
        
        if(magres_file->lattice == NULL) {
          return false;
        }

      } else if(strcmp(line->cols[0], "symmetry") == 0) {
        current_symmetry = magres_parse_symmetry(line);

        if(current_symmetry == NULL) {
          return false;
        }

        magres_file->symmetries[i_symmetry++] = current_symmetry;

      } else if(strcmp(line->cols[0], "units") == 0) {
        // Ignore units
      }
    }
  }

  return true;
}

MagresAtom *magres_find_atom(MagresFile *magres_file, const char* species, int index) {
  // Linear search for specific atom. It ain't pretty.

  int i;
  for(i=0; i<magres_file->num_atoms; ++i) {
    MagresAtom *atom = magres_file->atoms[i];

    if(strcmp(atom->species, species) == 0 && atom->index == index) {
      return atom;
    }
  }

  fprintf(stderr, "Could not find atom %s %d\n", species, index);

  return NULL;
}

MagresIsc *magres_parse_isc(MagresLine *line, MagresFile *magres_file) {
  if(line->num_cols != 14) {
    fprintf(stderr, "isc: Wrong number of columns, %d\n", line->num_cols);
    return NULL;
  }

  MagresIsc *isc = malloc(sizeof(MagresIsc));

  int species1_len = strlen(line->cols[1]);
  char *species1 = malloc(sizeof(char)*(species1_len + 1));
  memcpy(species1, line->cols[1], species1_len);
  species1[species1_len] = 0;
  
  int index1 = atoi(line->cols[2]);

  int species2_len = strlen(line->cols[3]);
  char *species2 = malloc(sizeof(char)*(species2_len + 1));
  memcpy(species2, line->cols[3], species2_len);
  species2[species2_len] = 0;
  
  int index2 = atoi(line->cols[4]);

  isc->atom1 = magres_find_atom(magres_file, species1, index1);
  isc->atom2 = magres_find_atom(magres_file, species2, index2);

  sscanf(line->cols[5], "%lg", &isc->K[0][0]);
  sscanf(line->cols[6], "%lg", &isc->K[0][1]);
  sscanf(line->cols[7], "%lg", &isc->K[0][2]);
  
  sscanf(line->cols[8], "%lg", &isc->K[1][0]);
  sscanf(line->cols[9], "%lg", &isc->K[1][1]);
  sscanf(line->cols[10], "%lg", &isc->K[1][2]);
  
  sscanf(line->cols[11], "%lg", &isc->K[2][0]);
  sscanf(line->cols[12], "%lg", &isc->K[2][1]);
  sscanf(line->cols[13], "%lg", &isc->K[2][2]);

  if(species1 != NULL) {
    free(species1);
    species1 = NULL;
  }

  if(species2 != NULL) {
    free(species2);
    species2 = NULL;
  }

  return isc;
}

MagresEfg *magres_parse_efg(MagresLine *line, MagresFile *magres_file) {
  if(line->num_cols != 12) {
    fprintf(stderr, "efg: Wrong number of columns, %d\n", line->num_cols);
    return NULL;
  }

  MagresEfg *efg = malloc(sizeof(MagresEfg));

  int species_len = strlen(line->cols[1]);
  char *species = malloc(sizeof(char)*(species_len + 1));
  memcpy(species, line->cols[1], species_len);
  species[species_len] = 0;
  
  int index = atoi(line->cols[2]);

  efg->atom = magres_find_atom(magres_file, species, index);

  sscanf(line->cols[3], "%lg", &efg->V[0][0]);
  sscanf(line->cols[4], "%lg", &efg->V[0][1]);
  sscanf(line->cols[5], "%lg", &efg->V[0][2]);
  
  sscanf(line->cols[6], "%lg", &efg->V[1][0]);
  sscanf(line->cols[7], "%lg", &efg->V[1][1]);
  sscanf(line->cols[8], "%lg", &efg->V[1][2]);
  
  sscanf(line->cols[9], "%lg", &efg->V[2][0]);
  sscanf(line->cols[10], "%lg", &efg->V[2][1]);
  sscanf(line->cols[11], "%lg", &efg->V[2][2]);

  if(species != NULL) {
    free(species);
    species = NULL;
  }

  return efg;
}

MagresMs *magres_parse_ms(MagresLine *line, MagresFile *magres_file) {
  if(line->num_cols != 12) {
    fprintf(stderr, "ms: Wrong number of columns, %d", line->num_cols);
    return NULL;
  }

  MagresMs *ms = malloc(sizeof(MagresMs));

  int species_len = strlen(line->cols[1]);
  char *species = malloc(sizeof(char)*(species_len + 1));
  memcpy(species, line->cols[1], species_len);
  species[species_len] = 0;
  
  int index = atoi(line->cols[2]);

  ms->atom = magres_find_atom(magres_file, species, index);

  sscanf(line->cols[3], "%lg", &ms->sigma[0][0]);
  sscanf(line->cols[4], "%lg", &ms->sigma[0][1]);
  sscanf(line->cols[5], "%lg", &ms->sigma[0][2]);
  
  sscanf(line->cols[6], "%lg", &ms->sigma[1][0]);
  sscanf(line->cols[7], "%lg", &ms->sigma[1][1]);
  sscanf(line->cols[8], "%lg", &ms->sigma[1][2]);
  
  sscanf(line->cols[9], "%lg", &ms->sigma[2][0]);
  sscanf(line->cols[10], "%lg", &ms->sigma[2][1]);
  sscanf(line->cols[11], "%lg", &ms->sigma[2][2]);

  if(species != NULL) {
    free(species);
    species = NULL;
  }

  return ms;
}


bool magres_parse_magres(MagresLine **lines, int num_lines, MagresFile *magres_file) {
  // Count how many of each line type we have and allocate space in the MagresFile
  int num_isc = 0, num_efg = 0, num_ms = 0;
  int i_line = 0;
  MagresLine *line = NULL;

  for(i_line = 0; i_line < num_lines; ++i_line) {
    line = lines[i_line];

    if(line->num_cols != 0) {
      if(strcmp(line->cols[0], "isc") == 0) {
        num_isc += 1;
      } else if(strcmp(line->cols[0], "efg") == 0) {
        num_efg += 1;
      } else if(strcmp(line->cols[0], "ms") == 0) {
        num_ms += 1;
      }
    }
  }

  if(num_isc > 0) {
    magres_file->isc = malloc(sizeof(MagresIsc*) * num_isc);
    magres_file->num_isc = num_isc;
  } else {
    magres_file->isc = NULL;
    magres_file->num_isc = 0;
  }

  if(num_efg > 0) {
    magres_file->efg = malloc(sizeof(MagresEfg*) * num_efg);
    magres_file->num_efg = num_efg;
  } else {
    magres_file->efg = NULL;
    magres_file->num_efg = 0;
  }

  if(num_ms > 0) {
    magres_file->ms = malloc(sizeof(MagresMs*) * num_ms);
    magres_file->num_ms = num_ms;
  } else {
    magres_file->ms = NULL;
    magres_file->num_ms = 0;
  }

  // Parse the lines into different types
  int i_isc=0, i_efg=0, i_ms=0;

  MagresIsc *current_isc = NULL;
  MagresEfg *current_efg = NULL;
  MagresMs *current_ms = NULL;

  for(i_line = 0; i_line < num_lines; ++i_line) {
    line = lines[i_line];

    if(line->num_cols != 0) {
      if(strcmp(line->cols[0], "isc") == 0) {
        current_isc = magres_parse_isc(line, magres_file);

        if(current_isc == NULL) {
          return false;
        }

        magres_file->isc[i_isc++] = current_isc;

      } else if(strcmp(line->cols[0], "efg") == 0) {
        current_efg = magres_parse_efg(line, magres_file);
        
        if(current_efg == NULL) {
          return false;
        }

        magres_file->efg[i_efg++] = current_efg;

      } else if(strcmp(line->cols[0], "ms") == 0) {
        current_ms = magres_parse_ms(line, magres_file);
        
        if(current_ms == NULL) {
          return false;
        }

        magres_file->ms[i_ms++] = current_ms;

      }
    }
  }

  return true;
}

int magres_parse_lines(char *block, MagresLine ***rtn_lines) {
  char *p = block;
  char *q = block;
  int line_size = 0;
  int num_lines = 0;

  // Count how many lines we have, for allocation
  while(true) {
    if(*p == '\n' || *p == 0) {
        line_size = p - q;

        if(line_size != 0) {
          num_lines += 1;
        }

        q = p + 1;
    }

    if(*p == 0) {
      break;
    } else {
     (++p);
    }
  }

  MagresLine **lines = NULL;
  MagresLine *line = NULL;

  lines = malloc(sizeof(MagresLine*) * num_lines);
  int i_line = 0;
  char *current_line = NULL;

  p = block;
  q = block;

  while(true) {
    if(*p == '\n' || *p == 0 || *p == '#') {
      if(current_line != NULL) {
        free(current_line);
        current_line = NULL;
      }

      line_size = p - q;

      if(line_size != 0) {
        current_line = malloc(sizeof(char) * (line_size + 1));
        memcpy(current_line, q, line_size);
        current_line[line_size] = 0;

        line = magres_tokenize(current_line);

        lines[i_line] = line;
        i_line++;
      }

      // If this is a comment, wind forward p to EOL or EOF
      if(*p == '#') {
        while(*p != '\n' && *p != 0) {
          ++p;
        }
      }

      q = p + 1;
    }

    if(*p == 0) {
      break;
    } else {
      (++p);
    }
  }
 
  if(current_line != NULL) {
    free(current_line);
    current_line = NULL;
  }

  *rtn_lines = lines;
  
  return num_lines;
}

// Parse a magres file from a string
bool magres_parse(MagresFile *magres_file, char *file)
{
  char *block_name = NULL;
  char *block_data = NULL;

  int block_name_size = 0;
  int block_data_size = 0;

  bool start_block = false;

  char *p = file;
 
  char *block_name_start = p;
  
  char *block_data_start = p;
  char *block_data_end = p;

  while(*p) {
    switch(*p) {
      case '#': // comment, skip rest of line
        while(*p != '\n' && *p != 0) {
          ++p;
        }

        break;

      case '[': // start of block tag
        if(*(p+1) == '/') { // end block tag
          if(!start_block) {
            fprintf(stderr, "ERROR: End block tag with no start block\n");
            return false;
          }

          start_block = false;
          block_data_end = p - 2;
        } else { // start block tag
          if(start_block) {
            fprintf(stderr, "ERROR: Start block tag inside block\n");
            return false;
          }

          block_name_start = p + 1;
          start_block = true;
        }
        break;

      case ']': // end of block tag
        if(start_block) {
          block_name_size = p - block_name_start;
          block_name = malloc(sizeof(char) * (block_name_size + 1));
          memcpy(block_name, block_name_start, block_name_size);
          block_name[block_name_size] = 0;

          //printf("Start of block %s\n", block_name);
          block_data_start = p + 1;
        } else {
          //printf("End of block %s\n", block_name);

          block_data_size = block_data_end - block_data_start + 1;
          block_data = malloc(sizeof(char) * (block_data_size + 1));
          memcpy(block_data, block_data_start, block_data_size);
          block_data[block_data_size] = 0;

          int num_lines = 0;
          MagresLine **lines = NULL;

          num_lines = magres_parse_lines(block_data, &lines);

          if(strcmp(block_name, "atoms") == 0) {
            bool rtn = magres_parse_atoms(lines, num_lines, magres_file);
            
            if(!rtn) {
              return false;
            }
          } else if(strcmp(block_name, "magres") == 0) {
            bool rtn = magres_parse_magres(lines, num_lines, magres_file);

            if(!rtn) {
              return false;
            }
          }
  
          // Deallocation
          int i_line = 0;
          for(; i_line < num_lines; ++i_line) {
            if(lines[i_line] != NULL) {
              magres_line_dealloc(lines[i_line]);
              free(lines[i_line]);
              lines[i_line] = NULL;
            }
          }
          
          if(lines != NULL) {
            free(lines);
            lines = NULL;
          }

          free(block_name);
          block_name = NULL;
        }
        break;
    }

    (++p);
  }

  if(block_name != NULL) {
    printf("Unterminated block %s\n", block_name);
    free(block_name);
    block_name = NULL;
  }

  return true;
}

char *read_file(const char* path) {
  FILE* fp = fopen(path, "r");

  char *buffer;

  if(fp != NULL) {
    if(fseek(fp, 0L, SEEK_END) == 0) {
      int buffer_size = ftell(fp);

      if(buffer_size == -1) return NULL;

      buffer = malloc(sizeof(char) * (buffer_size + 1));

      if(fseek(fp, 0L, SEEK_SET) != 0) return NULL;

      size_t new_length = fread(buffer, sizeof(char), buffer_size, fp);
      if(new_length == 0) {
        fputs("Error reading file\n", stderr);
      } else {
        buffer[++new_length] = '\0'; /* Just to be safe. */
      }
    }
  } else {
    return NULL;
  }

  return buffer;
}

int main(int argc, const char **argv) {
  char *magres_file_data = read_file(argv[1]);
  
  MagresFile *magres_file = malloc(sizeof(MagresFile));

  magres_file_init(magres_file);

  if(magres_file_data != NULL ) {
    bool rtn = magres_parse(magres_file, magres_file_data);

    if(!rtn) {
      printf("Error parsing file\n");
      return 0;
    }

  } else {
    printf("Error loading file\n");
  }

  if(magres_file_data != NULL) {
    free(magres_file_data);
    magres_file_data = NULL;
  }

  printf("%d atoms\n", magres_file->num_atoms);
  printf("%d symmetries\n", magres_file->num_symmetries);

  if(magres_file->lattice != NULL) {
    printf("Has lattice\n");
  } else {
    printf("No lattice\n");
  }

  printf("%d J-couplings\n", magres_file->num_isc);
  printf("%d EFG tensors\n", magres_file->num_efg);
  printf("%d MS tensors\n", magres_file->num_ms);

  printf("Atoms:\n");
  int i;
  MagresAtom *atom = NULL;
  for(i = 0; i<magres_file->num_atoms; ++i) {
    atom = magres_file->atoms[i];
    printf("  %i %s %s %f %f %f\n", atom->index, atom->species, atom->label, atom->position[0],  atom->position[1], atom->position[2]);
  }

  printf("Lattice:\n");
  printf("  %f %f %f\n", magres_file->lattice->lattice[0][0],  magres_file->lattice->lattice[0][1], magres_file->lattice->lattice[0][2]);
  printf("  %f %f %f\n", magres_file->lattice->lattice[1][0],  magres_file->lattice->lattice[1][1], magres_file->lattice->lattice[1][2]);
  printf("  %f %f %f\n", magres_file->lattice->lattice[2][0],  magres_file->lattice->lattice[2][1], magres_file->lattice->lattice[2][2]);

  printf("Symmetries:\n");
  for(i = 0; i<magres_file->num_symmetries; ++i) {
    printf("%s\n", magres_file->symmetries[i]->symmetry_string);
  }

  MagresIsc *isc = NULL;
  for(i = 0; i<magres_file->num_isc; ++i) {
    isc = magres_file->isc[i];

    double K_iso = (isc->K[0][0] + isc->K[1][1] + isc->K[2][2])/3.0;
    printf("ISC: %s %d --> %s %d = %f\n", isc->atom1->species, isc->atom1->index, isc->atom2->species, isc->atom2->index, K_iso);
  }

  MagresMs *ms = NULL;
  for(i = 0; i<magres_file->num_ms; ++i) {
    ms = magres_file->ms[i];

    double ms_iso = (ms->sigma[0][0] + ms->sigma[1][1] + ms->sigma[2][2])/3.0;
    printf("MS: %s %d = %f\n", ms->atom->species, ms->atom->index, ms_iso);
  }

  if(magres_file != NULL) {
    magres_file_dealloc(magres_file);
    free(magres_file);
    magres_file = NULL;
  }
}
