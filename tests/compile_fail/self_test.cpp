// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Checks the compile-fail test mechanism itself: this file compiles as
// written and fails when CANARY_COMPILE_FAIL is defined.
#ifdef CANARY_COMPILE_FAIL
static_assert(false, "compile-fail self test");
#endif
