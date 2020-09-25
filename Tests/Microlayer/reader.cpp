template <class T>
bool read_variable(std::ifstream & input, T & var) {
  std::string line;
  if(getline(input,line)) {
    std::istringstream ssline(line);
    ssline >> var;
    OPR(var);
    return 0;
  } else {
    return 1;
  }
}
