#ifndef CANARY_DATABASE_FILE_IO_HPP
#define CANARY_DATABASE_FILE_IO_HPP

#include <filesystem>
#include <istream>
#include <string>
#include <string_view>

namespace canary::database::detail
{
// Reads through EOF; a read failure discards partial contents and reports the path.
std::string
  read_contents(std::istream& input, const std::filesystem::path& path, std::string_view kind);
std::string read_text(const std::filesystem::path& path, std::string_view kind);

// Commits a fully serialized document by renaming a completed neighbouring file.
// Rejects symlink/non-regular destinations and preserves existing permission bits.
void write_text(const std::filesystem::path& path, const std::string& text);
}  // namespace canary::database::detail

#endif  // CANARY_DATABASE_FILE_IO_HPP
