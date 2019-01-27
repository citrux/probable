#pragma once
#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_PROBABLE
    #define PROBABLE_PUBLIC __declspec(dllexport)
  #else
    #define PROBABLE_PUBLIC __declspec(dllimport)
  #endif
#else
  #ifdef BUILDING_PROBABLE
      #define PROBABLE_PUBLIC __attribute__ ((visibility ("default")))
  #else
      #define PROBABLE_PUBLIC
  #endif
#endif

namespace probable {

class PROBABLE_PUBLIC Probable {

public:
  Probable();
  int get_number() const;

private:

  int number;

};

}

