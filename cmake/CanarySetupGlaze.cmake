include(CPM)

# JSON reading and writing for the species database library only.
CPMAddPackage(
    NAME glaze
    GITHUB_REPOSITORY stephenberry/glaze
    GIT_TAG v9.0.0
    GIT_SHALLOW YES
    SYSTEM YES
    EXCLUDE_FROM_ALL YES
    OPTIONS
        "glaze_DEVELOPER_MODE OFF"
        "glaze_BUILD_EXAMPLES OFF"
        "glaze_ENABLE_FUZZING OFF"
)
