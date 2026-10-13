// Writes the JSON Schema of collision-integral files to standard output.

#include <cstdio>
#include <canary/database.hpp>

int main()
{
    const auto schema = canary::database::collision_integrals_schema();
    std::fputs(schema.c_str(), stdout);
}
