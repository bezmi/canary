// Writes the JSON Schema of species files to standard output.

#include <cstdio>
#include <gaslib/database.hpp>

int main()
{
    const auto schema = gaslib::database::species_schema();
    std::fputs(schema.c_str(), stdout);
}
