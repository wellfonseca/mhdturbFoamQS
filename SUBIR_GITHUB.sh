#!/bin/bash
# SUBIR_GITHUB.sh <utilizador-github> [nome-do-repo]
# Configura o remote e envia o repositorio para o GitHub.
set -e
U=${1:?uso: bash SUBIR_GITHUB.sh <utilizador-github> [nome-do-repo]}
R=${2:-mhdturbFoamQS}
cd "$(dirname "$0")" || exit 1
echo "=== a testar a autenticacao ==="
ssh -T git@github.com 2>&1 | head -2 || true
git remote remove origin 2>/dev/null || true
git remote add origin "git@github.com:$U/$R.git"
echo "=== a enviar para git@github.com:$U/$R.git ==="
git push -u origin main
echo
echo "feito. Repositorio em https://github.com/$U/$R"
