// Writes the JSON Schema of reactions files to standard output.

#include <cstdio>
#include <canary/database.hpp>

int main()
{
    const auto schema = canary::database::reactions_schema();
    std::fputs(schema.c_str(), stdout);
}
