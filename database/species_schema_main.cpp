// Writes the JSON Schema of species files to standard output.

#include <cstdio>
#include <canary/database.hpp>

int main()
{
    const auto schema = canary::database::species_schema();
    std::fputs(schema.c_str(), stdout);
}
