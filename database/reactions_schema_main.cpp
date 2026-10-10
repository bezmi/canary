// Writes the JSON Schema of reactions files to standard output.

#include <cstdio>
#include <gaslib/database.hpp>

int main()
{
    const auto schema = gaslib::database::reactions_schema();
    std::fputs(schema.c_str(), stdout);
}
