#include <iostream>
#include <probable.hh>

int main(int argc, char **argv) {
  if (argc != 1) {
    std::cout << argv[0] << " takes no arguments.\n";
    return 1;
  }
  double c = probable::uniform();
  std::cout << c << "\n";
  return !(c < 1);
}
