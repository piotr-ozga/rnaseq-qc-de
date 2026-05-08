#!/usr/bin/env bash
set -euo pipefail

TESTS=()

# Always run the pipeline level integration test
[[ -f "tests/main.nf.test" ]] && TESTS+=("tests/main.nf.test")

for FILE in ${CHANGED_FILES:-}; do
  # Test file modified directly
  if [[ "$FILE" == *.nf.test && -f "$FILE" ]]; then
    TESTS+=("$FILE")
    
  # Local module modification
  elif [[ "$FILE" == modules/local/* ]]; then
    MODULE_PATH=$(echo "$FILE" | grep -oP 'modules/local/[^/]+')
    while IFS= read -r TEST; do
      TESTS+=("$TEST")
    done < <(find "tests/${MODULE_PATH}" -name "main.nf.test" 2>/dev/null)

  # Local subworkflow modification
  elif [[ "$FILE" == subworkflows/local/* ]]; then
    SUBWORKFLOW=$(echo "$FILE" | cut -d'/' -f3)
    TEST="tests/subworkflows/local/${SUBWORKFLOW}/main.nf.test"
    if [[ -f "$TEST" ]]; then
      TESTS+=("$TEST")
    fi
  fi
done

# Deduplicate and convert to JSON array for GitHub Actions matrix
MATRIX=$(printf '%s\n' "${TESTS[@]:-}" | sort -u | jq -R . | jq -sc .)
echo "matrix=$MATRIX" >> "$GITHUB_OUTPUT"