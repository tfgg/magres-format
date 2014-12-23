
#ifdef USE_LIBCMATRIX
#include "libcmatrix/geometry.h"
using namespace libcmatrix;
typedef vector3 vector_t;
typedef rmatrix3 tensor_t;
#include "libcmatrix/List.h"
#define LIST_T libcmatrix::List
typedef Failed exception_t;
typedef InvalidParameter notmagres_exception_t;
#else
typedef double vector_t [3];
typedef double tensor_t [3][3];
#include <list>
#define LIST_T std::list
typedef std::runtime_error exception_t;
typedef std::runtime_error notmagres_exception_t; //!< should use different type so can distinguish failure modes
#endif

// size of buffer to allocate for error messages
#define ERRM_MAXLEN 1024

static const char line_separator='\n';
static const char token_separator[] = " \t\r";
#define MAJOR_VERSION_LIMIT 1

#include <string>
#include <string.h>

// The atomic properties

struct atomid_t {
  atomid_t(const char* speciesv, size_t indexv)
    : species(speciesv), index(indexv) {}

  std::string species;
  size_t index;
};

struct MagresLine;

struct MagresAtom {
  MagresAtom(const MagresLine&);
  size_t index;
  std::string species; 
  std::string label;
  vector_t position;

  std::ostream& printname(std::ostream& ostr) const { return ostr << species << ' ' << index; };
  friend std::ostream& operator<< (std::ostream& ostr, const MagresAtom& atom);

  bool operator== (const atomid_t& id) const 
  { return (id.index==index) && (strcmp(species.c_str(),id.species.c_str())==0); }

};

struct MagresLattice {
  tensor_t lattice;
  MagresLattice(const MagresLine&);
};

struct MagresSymmetry {
  std::string symmetry_string;
  MagresSymmetry(const MagresLine&);
};

struct MagresFile;

// The magnetic resonance properties
struct MagresMs {
  const MagresAtom *atom;
  tensor_t sigma;
  MagresMs(const MagresLine&, MagresFile&);
};

struct MagresIsc {
  const MagresAtom* atom1;
  const MagresAtom* atom2;
  tensor_t K;
  MagresIsc(const MagresLine&, MagresFile&);
};

struct MagresEfg {
  const MagresAtom* atom;
  tensor_t V;
  MagresEfg(const MagresLine&, MagresFile&);
};

enum block_t { ATOM, MAGRES };

// The whole file
struct MagresFile {
  MagresFile() : latticep(NULL) {}
  ~MagresFile() { delete latticep; }

  MagresLattice* latticep;
  
  size_t num_atoms() const { return atoms.size(); }
  LIST_T<MagresAtom> atoms;

  size_t num_symmetries() const { return symmetries.size(); }
  LIST_T<MagresSymmetry> symmetries;

  size_t num_ms() const { return ms.size(); }
  LIST_T<MagresMs> ms;

  size_t num_isc() const { return isc.size(); }
  LIST_T<MagresIsc> isc;

  size_t num_efg() const { return efg.size(); }
  LIST_T<MagresEfg> efg;

  const MagresAtom* find_atom(const atomid_t&) const;
  const MagresAtom* parse_find_atom(const MagresLine& line,size_t base) const;
  void parse_from_string(char* file);
  void parse_from_file(const char* file);
  void parse_atom(const MagresLine&);
  void parse_magres(const MagresLine&);
  void parse_lines(char *block, block_t blocktype);
  char* windforward(char*);

  friend std::ostream& operator<< (std::ostream& ostr, const MagresFile& magres_file);
};

// Generic line tokenizer
struct MagresLine {
  MagresLine(char*);

  size_t num_cols() const { return cols.size(); }
  const char* column(size_t i) const { return cols(i).c_str(); }
  void verify_columns(size_t cols, const char* label) const;
  void read_tensor(tensor_t& dest, size_t base) const;

  LIST_T<std::string> cols;
};

MagresLine::MagresLine(char* p)
{
  char* token;
  while ((token=strtok(p,token_separator))) {
    p=NULL;
    cols.push_back(token);
  }
}

long parse_integer(const char* p)
{
  char* tail;
  long val=strtol(p,&tail,10);
  if (*tail!='\0') {
    char errm[ERRM_MAXLEN];
    sprintf(errm,"failed to parse %s as integer", p);
    throw exception_t(errm);
  }
  return val;
}

double parse_double(const char* p)
{
  char* tail;
  double val=strtod(p,&tail);
  if (*tail!='\0') {
    char errm[ERRM_MAXLEN];
    sprintf(errm,"failed to parse %s as floating point value", p);
    throw exception_t(errm);
  }
  return val;
}

void MagresLine::verify_columns(size_t cols, const char* label) const
{
  if (num_cols() != cols) {
    char errm[ERRM_MAXLEN];
    sprintf(errm,"%s: Wrong number of columns, %d", label, num_cols());
    throw exception_t(errm);
  }
}

// [atoms] block parser
MagresAtom::MagresAtom(const MagresLine& line)
{
  line.verify_columns(7,"atom");
  species=line.column(1);
  label=line.column(2);
  index = parse_integer(line.column(3));
  position.x = parse_double(line.column(4));
  position.y = parse_double(line.column(5));
  position.z = parse_double(line.column(6));
}

MagresSymmetry::MagresSymmetry(const MagresLine& line)
{
  line.verify_columns(2,"symmetry");
  symmetry_string=line.column(1);
}

void MagresLine::read_tensor(tensor_t& dest, size_t base) const
{
  dest[0][0] = parse_double(column(base));
  dest[0][1] = parse_double(column(base+1));
  dest[0][2] = parse_double(column(base+2));

  dest[1][0] = parse_double(column(base+3));
  dest[1][1] = parse_double(column(base+4));
  dest[1][2] = parse_double(column(base+5));

  dest[2][0] = parse_double(column(base+6));
  dest[2][1] = parse_double(column(base+7));
  dest[2][2] = parse_double(column(base+8));  
}

MagresLattice::MagresLattice(const MagresLine& line)
{
  line.verify_columns(10,"lattice");
  line.read_tensor(lattice,1);
}

void MagresFile::parse_atom(const MagresLine& line)
{
  const char* id=line.column(0);
  if (strcmp(id, "atom") == 0)
    atoms.push_back(MagresAtom(line));
  else if (strcmp(id, "lattice") == 0)
    latticep= new MagresLattice(line);
  else if (strcmp(id, "symmetry") == 0)
    symmetries.push_back(line);
  //ignore everything else, including units
}

const MagresAtom* MagresFile::find_atom(const atomid_t& searchatom) const
{
  LIST_T<MagresAtom>::const_iterator iter=std::find(atoms.begin(),atoms.end(),searchatom);
  if (iter!=atoms.end())
    return &(*iter);

  fprintf(stderr, "Could not find atom %s %d", searchatom.species.c_str(), searchatom.index);
  return NULL;
}

const MagresAtom* MagresFile::parse_find_atom(const MagresLine& line,size_t base) const
{
  const atomid_t searchatom( line.column(base), parse_integer(line.column(base+1)));
  return find_atom(searchatom);
}

MagresIsc::MagresIsc(const MagresLine& line, MagresFile& magres_file) 
{
  line.verify_columns(14,"isc");

  atom1=magres_file.parse_find_atom(line,1);
  atom2=magres_file.parse_find_atom(line,3);

  line.read_tensor(K,5);
}

MagresEfg::MagresEfg(const MagresLine& line, MagresFile& magres_file)
{
  line.verify_columns(12,"efg");
  atom = magres_file.parse_find_atom(line,1);
  line.read_tensor(V,3);
}

MagresMs::MagresMs(const MagresLine& line, MagresFile& magres_file)
{
  line.verify_columns(12,"ms");
  atom = magres_file.parse_find_atom(line,1);
  line.read_tensor(sigma,3);
}

void MagresFile::parse_magres(const MagresLine& line)
{
  const char* id=line.column(0);
  
  if (strcmp(id, "isc") == 0)
    isc.push_back(MagresIsc(line,*this));
  else if (strcmp(id, "efg") == 0)
    efg.push_back(MagresEfg(line,*this));
  else if (strcmp(id, "ms") == 0)
    ms.push_back(MagresMs(line,*this));
}

char* MagresFile::windforward(char* p)
{
  while( (*p != line_separator) && (*p != 0))
    ++p;
  return p;
}

void MagresFile::parse_lines(char *block, block_t blocktype)
{
  char* lastp=block;

  for (char* p=block;;p++) {
    if ( (*p == line_separator) || (*p == '\0') || (*p == '#') ) {
      if (p>lastp) {
	char storep=*p;	
	*p='\0'; // terminate
        MagresLine line(lastp);
	if (line.num_cols()) {
	  switch (blocktype) {
	  case ATOM:
	    parse_atom(line);
	    break;
	  case MAGRES:
	    parse_magres(line);
	    break;
	  }
	}
	*p=storep; //!< restore
      }

      // If this is a comment, wind forward p to EOL or EOF
      if (*p == '#')
	p=windforward(p);

      lastp = p + 1;
      if (*p == '\0') 
	break;
    }
  }
}

// Parse a magres file from a string
void MagresFile::parse_from_string(char* file)
{
  bool start_block = false; 
  char* block_name = NULL;
  char* block_data_end = NULL;
  char* block_data_start = NULL;

  unsigned int major,minor;
  if (sscanf(file,"#$magres-abinitio-v%u.%u",&major,&minor)!=2)
    throw notmagres_exception_t("ERROR: Not a new-style magres file");
  
  if (major>MAJOR_VERSION_LIMIT) {
    char errm[ERRM_MAXLEN];
    sprintf(errm,"ERROR: magres version (%d.%d) exceeds limit (%d)",major,minor,MAJOR_VERSION_LIMIT);
    throw exception_t(errm);
  }

  for (char* p= file;*p;p++) {
    switch (*p) {
      case '#': // comment, skip rest of line
        p=windforward(p);
        break;
	
    case '[': // start of block tag
      if(*(p+1) == '/') { // end block tag
	if(!start_block)
	  throw exception_t("ERROR: End block tag with no start block");
	
	start_block = false;
	block_data_end = p - 2;
      } else { // start block tag
	if (start_block)
	  throw exception_t("ERROR: Start block tag inside block");
	
	block_name = p + 1;
	start_block = true;
      }
      break;
      
    case ']': // end of block tag
      if (start_block) {
	*p='\0';
	block_data_start = p + 1;
      } 
      else {
	if ((block_data_end==NULL) || (block_data_start==NULL))
	  throw exception_t("ERROR: block data start/end unset");
	block_data_end[1]='\0'; //!< terminate data block
	if (strcmp(block_name, "atoms") == 0)
	  parse_lines(block_data_start,ATOM);
	else if (strcmp(block_name, "magres") == 0)
	  parse_lines(block_data_start,MAGRES);
	block_name=NULL;
      }
    }
  }
  
  if (block_name != NULL) {
    char errm[ERRM_MAXLEN];
    sprintf(errm,"Unterminated block %s", block_name);
    throw exception_t(errm);
  }
}

 class FileOpenGuard {
 public:
   FileOpenGuard(const char*);
   ~FileOpenGuard() { fclose(fp); }
   FILE* operator()() { return fp; }
 private:
   FILE* fp;
 };

 FileOpenGuard::FileOpenGuard(const char* path)
   {
     fp = fopen(path, "r");          
     char errm[ERRM_MAXLEN];
     if (fp==NULL) {
       sprintf(errm,"Failed to open file for reading: %s",path);
       throw exception_t(errm);
     }
   }
 

 void MagresFile::parse_from_file(const char* path)
   {
     FileOpenGuard FP(path);
     FILE* fp=FP();
     std::string contents;     
     if (fseek(fp, 0L, SEEK_END) == 0) {
       const int buffer_size = ftell(fp);       
       if (buffer_size == -1)
	 throw exception_t("Error reading file size");
       
       contents.resize(buffer_size+1);
       std::rewind(fp);
       const size_t new_length = std::fread(&contents[0], 1, buffer_size, fp);
       if(new_length == 0)
	 throw exception_t("Error reading file");
       
       contents[new_length+1] = '\0'; /* Just to be safe. */
       char* asraw=const_cast<char*>(contents.c_str()); //marginally dodgy - need to assume that parse_from_string doesn't over-run
       parse_from_string(asraw);
     }
     else
       throw exception_t("Error finding file end");
   }

std::ostream& print_tensor(const tensor_t& A, std::ostream& ostr)
   {     
     ostr << "  " << A[0][0] << ' ' << A[0][1] << ' ' << A[0][2] << '\n';
     ostr << "  " << A[1][0] << ' ' << A[1][1] << ' ' << A[1][2] << '\n';
     ostr << "  " << A[2][0] << ' ' << A[2][1] << ' ' << A[2][2] << '\n';
     return ostr;
   }
 
 std::ostream& operator<< (std::ostream& ostr, const MagresAtom& atom)
   {
     return ostr << "  " << atom.index << ' ' << atom.species <<  ' ' << atom.label <<  ' ' << atom.position << '\n';
   }

 std::ostream& operator<< (std::ostream& ostr, const MagresFile& magres_file)
   {
     ostr << magres_file.num_atoms() << " atoms\n";
     ostr << magres_file.num_symmetries() << " symmetries\n";
     
     if (magres_file.latticep) {
       ostr <<  "Lattice:\n";
       print_tensor((magres_file.latticep)->lattice,ostr);
     }
     else
       ostr << "No lattice\n";

     ostr << magres_file.num_isc() << " J-couplings\n";
     ostr << magres_file.num_efg() << " EFG tensors\n";
     ostr << magres_file.num_ms() << " MS tensors\n";

     ostr << "Atoms:\n";
     for (size_t i = 0; i<magres_file.num_atoms(); ++i)
       ostr << magres_file.atoms(i);
     
     printf("Symmetries:\n");
     for (size_t i = 0; i<magres_file.num_symmetries(); ++i)
       ostr << magres_file.symmetries(i).symmetry_string << '\n';

     for (size_t i = 0; i<magres_file.num_isc(); ++i) {
       const MagresIsc& isc(magres_file.isc(i));
       const double K_iso = (isc.K[0][0] + isc.K[1][1] + isc.K[2][2])/3.0;
       ostr << "ISC: ";
       (isc.atom1)->printname(ostr) << " --> ";
       (isc.atom2)->printname(ostr) << " = " << K_iso << '\n';
     }

     for (size_t i = 0; i<magres_file.num_ms(); ++i) {
       const MagresMs& ms(magres_file.ms(i));       
       const double ms_iso = (ms.sigma[0][0] + ms.sigma[1][1] + ms.sigma[2][2])/3.0;
       ostr << "MS: ";
       (ms.atom)->printname(ostr) << " = " << ms_iso << '\n';
     }

     return ostr;
   }

int main(int argc, const char **argv) {
  try {
    MagresFile magres;
    magres.parse_from_file(argv[1]);
    std::cout << magres;
  }
  catch (exception_t& exc) {
    std::cerr << "Parsing failed: " << exc << '\n';
    return 1;
  }
  catch (notmagres_exception_t&) {
    std::cerr << "Not a new format magres file\n";
    return 2;
  }
  return 0;
}
