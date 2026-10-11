// Writes the JSON Schema of energy-exchange files to standard output.

#include <cstdio>
#include <gaslib/database.hpp>

int main()
{
    const auto schema = gaslib::database::energy_exchange_schema();
    std::fputs(schema.c_str(), stdout);
}
