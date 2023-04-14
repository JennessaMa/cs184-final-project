#include "CGL/CGL.h"
#include "CGL/viewer.h"
// #include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#define TINYEXR_IMPLEMENTATION
#include "CGL/tinyexr.h"
// typedef uint32_t gid_t;
#include <iostream>

#include "svg.h"
#include "drawrend.h"
#include "transforms.h"
#include "svgparser.h"

#include <ft2build.h>
#include FT_FREETYPE_H


using namespace std;
using namespace CGL;

#define msg(s) cerr << "[Drawer] " << s << endl;


SVG *loadFile( const char* path ) {

  SVG* svg = new SVG();

  if( SVGParser::load( path, svg ) < 0) {
    delete svg;
    return NULL;
  }
  
  return svg;
}

vector<SVG*> loadDirectory( const char* path ) {

  DIR *dir = opendir (path);
  vector<SVG*> svgs;
  if(dir) {
    
    struct dirent *ent; size_t n = 0;
    
    // load files
    string pathname = path; 
    if (pathname[pathname.size()-1] != '/') pathname.push_back('/');
    while (((ent = readdir (dir)) != NULL) && (n < 9)) {

      string filename = ent->d_name;
      string filesufx = filename.substr(filename.find_last_of(".") + 1, 3);
      if (filesufx == "svg" ) {
        cerr << "[Drawer] Loading " << filename << "... "; 
        SVG *svg = loadFile((pathname + filename).c_str());
        if (!svg) {
          cerr << "Failed (Invalid SVG file)" << endl;
        } else {
          cerr << "Succeeded" << endl;
          svgs.push_back(svg);
          n++;
        }
      }
    }

    closedir (dir);

    if (n) {
      msg("Successfully Loaded " << n << " files from " << path);
      return svgs;
    }

    msg("No valid svg files found in " << path);
    return svgs;
  } 

  msg("Could not open directory" << path);
  return svgs;
}

vector<SVG*> loadPath( const char* path ) {

  struct stat st;

  // file exist?
  if(stat(path, &st) < 0 ) {
    msg("File does not exit: " << path);
    return vector<SVG*>();
  }

  // load directory
  if( st.st_mode & S_IFDIR ) {
    return loadDirectory(path);
  } 

  // load file
  if( st.st_mode & S_IFREG ) {
    return vector<SVG*>(1,loadFile(path));
  }

  msg("Invalid path: " << path);
  return vector<SVG*>();
}


int main( int argc, char** argv ) {

  if (argc < 2) {
    msg("Not enough arguments. Pass in an .svg or a directory of .svg files.");
    return 0;
  }

  vector<SVG*> svgs(loadPath(argv[1]));
  if (svgs.empty()) {
    msg("No svg files successfully loaded. Exiting.");
    return 0;
  }

  FT_Library ft;
  FT_Face face;
  if (FT_Init_FreeType (&ft) == 0) {
    cout << "Successfully initialized freetype :)" << endl;
    auto error = FT_New_Face( ft,
                         "../15716_FrugalSansLight.ttf",
                         0,
                         &face );
    if (error == FT_Err_Unknown_File_Format) {
      cout << "failed to open ttf" << endl;
    } else if (error) {
      cout << "some other error while opening ttf" << endl;
      cout << FT_Error_String(error) << endl;
    }
    FT_Set_Char_Size(face, 0, 16 * 64, 300, 300);
    error = FT_Load_Glyph(face, FT_Get_Char_Index(face, 'a'), FT_LOAD_DEFAULT);
    FT_Outline *outline = &face->glyph->outline;

    cout << outline->n_points << endl;
    for (int i = 0; i < outline->n_points; i += 1) {
      cout << outline->points[i].x << endl;
      cout << outline->points[i].y << endl;
      // rasterize_point(outline->points[i].x, outline->points[i].y, Color(1, 1, 1));
    }

    // FT_Render_Glyph(face->glyph, ft_render_mode_normal);
  } else {
    cout << "Failed to initialize freetype :(" << endl;
  }

//
//    // create application
//  DrawRend app(svgs);
//
//  if (argc > 4 && strcmp(argv[2],"nogl") == 0) {
//    app.init();
//    app.set_gl(false);
//    app.resize(stoi(argv[3]), stoi(argv[4]));
//    app.write_framebuffer();
//    return 0;
//  }
//
//  // create viewer
//  Viewer viewer = Viewer();
//
//  // set renderer
//  viewer.set_renderer(&app);
//
//  // init viewer
//  viewer.init();
//
//  // start viewer
//  viewer.start();

  exit(EXIT_SUCCESS);

  return 0;

}
