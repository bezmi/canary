// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#include "file_io.hpp"

#include <algorithm>
#include <array>
#include <canary/database.hpp>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <optional>
#include <random>
#include <string>
#include <system_error>
#include <utility>

namespace canary::database::detail
{
namespace
{
std::optional<std::filesystem::perms> destination_permissions(const std::filesystem::path& path)
{
    std::error_code error;
    const auto status = std::filesystem::symlink_status(path, error);
    if (error && error != std::errc::no_such_file_or_directory) {
        throw DatabaseError(
          std::format("cannot inspect output file {}: {}", path.string(), error.message()));
    }
    if (status.type() == std::filesystem::file_type::not_found) { return {}; }
    if (!std::filesystem::is_regular_file(status)) {
        throw DatabaseError(
          std::format("output file {} must be a regular file, not a symlink or other file type",
                      path.string()));
    }
    return status.permissions();
}

void preserve_permissions(const std::filesystem::path& temporary,
                          const std::filesystem::path& destination,
                          std::optional<std::filesystem::perms> permissions)
{
    if (!permissions) { return; }
    std::error_code error;
    std::filesystem::permissions(temporary,
                                 *permissions,
                                 std::filesystem::perm_options::replace,
                                 error);
    if (error) {
        throw DatabaseError(std::format("cannot preserve permissions of output file {}: {}",
                                        destination.string(),
                                        error.message()));
    }
}

// Owns only a file successfully created with exclusive access. Cleanup runs on
// every failure before replacement; committed files are released from ownership.
class TemporaryOutput
{
  public:
    explicit TemporaryOutput(const std::filesystem::path& destination)
    {
        auto parent = destination.parent_path();
        if (parent.empty()) { parent = "."; }
        std::random_device random;
        for (int attempt = 0; attempt < 128; ++attempt) {
            auto candidate = parent / (".canary-" + std::to_string(random()) + "-" +
                                       std::to_string(random()) + ".tmp");
            stream.open(candidate, std::ios::binary | std::ios::out | std::ios::noreplace);
            if (stream.is_open()) {
                path = std::move(candidate);
                return;
            }
            std::error_code error;
            const auto status = std::filesystem::symlink_status(candidate, error);
            if (error || !std::filesystem::exists(status)) {
                throw DatabaseError(
                  std::format("cannot create temporary output for {}", destination.string()));
            }
            stream.clear();
        }
        throw DatabaseError(
          std::format("cannot create unique temporary output for {}", destination.string()));
    }
    TemporaryOutput(const TemporaryOutput&) = delete;
    TemporaryOutput& operator=(const TemporaryOutput&) = delete;

    ~TemporaryOutput()
    {
        if (stream.is_open()) { stream.close(); }
        if (!path.empty()) {
            std::error_code error;
            std::filesystem::remove(path, error);
        }
    }

    std::filesystem::path path;
    std::ofstream stream;
};
}  // namespace

std::string
  read_contents(std::istream& input, const std::filesystem::path& path, std::string_view kind)
{
    std::string text;
    std::array<char, 16384> buffer;
    for (;;) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        if (input.bad() || (input.fail() && !input.eof())) {
            throw DatabaseError(std::format("cannot read {} file {}", kind, path.string()));
        }
        text.append(buffer.data(), static_cast<std::size_t>(input.gcount()));
        if (input.eof()) { return text; }
    }
}

std::string read_text(const std::filesystem::path& path, std::string_view kind)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) { throw DatabaseError(std::format("cannot open {} file {}", kind, path.string())); }
    return read_contents(input, path, kind);
}

void write_text(const std::filesystem::path& path, const std::string& text)
{
    const auto initial_permissions = destination_permissions(path);
    TemporaryOutput pending(path);
    // Restrict access before writing any contents; reapply after close because
    // writing can clear special permission bits on some filesystems.
    preserve_permissions(pending.path, path, initial_permissions);
    // Bound each stream operation without narrowing the document's total size.
    for (std::size_t offset = 0; offset < text.size();) {
        const auto count = std::min<std::size_t>(text.size() - offset, 16384);
        pending.stream.write(text.data() + offset, static_cast<std::streamsize>(count));
        if (!pending.stream) {
            throw DatabaseError(std::format("cannot write output file {}", path.string()));
        }
        offset += count;
    }
    pending.stream.flush();
    if (!pending.stream) {
        throw DatabaseError(std::format("cannot flush output file {}", path.string()));
    }
    pending.stream.close();
    if (!pending.stream) {
        throw DatabaseError(std::format("cannot close output file {}", path.string()));
    }

    // Recheck after writing so existing destination permissions are copied as
    // late as possible. Renaming replaces the directory entry, not a link target.
    preserve_permissions(pending.path, path, destination_permissions(path));
    std::error_code error;
    std::filesystem::rename(pending.path, path, error);
    if (error) {
        throw DatabaseError(
          std::format("cannot replace output file {}: {}", path.string(), error.message()));
    }
    pending.path.clear();
}
}  // namespace canary::database::detail
