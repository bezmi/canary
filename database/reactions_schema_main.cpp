// Writes the JSON Schema of reactions files to standard output.

#include <cstdio>
#include <canary/database.hpp>

int main()
{
    try {
        const auto schema = canary::database::reactions_schema();
        if (std::fputs(schema.c_str(), stdout) == EOF || std::fflush(stdout) == EOF) {
            std::fputs("cannot write schema to standard output\n", stderr);
            return 1;
        }
    } catch (const canary::database::DatabaseError& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
