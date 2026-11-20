// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Writes the JSON Schema of reactions files to standard output.

#include <canary/database.hpp>
#include <cstdio>

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
