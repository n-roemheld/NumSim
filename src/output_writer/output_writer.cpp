#include "output_writer/output_writer.h"

#include <iostream>
#include <string>
#include <cstring>

OutputWriter::OutputWriter(std::shared_ptr<Discretization> discretization, std::string outputFolder)
 : discretization_(discretization), fileNo_(0), outputFolder_(outputFolder)
{
  // create "out" subdirectory if it does not yet exist
  std::string commandString = "mkdir -p" + outputFolder_;
  char * commandChar = new char [commandString.length() + 1];
  strcpy(commandChar, commandString.c_str());
  int returnValue = system(commandChar); // Hilfe
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \" << outputFolder_ << \"." << std::endl;
}
