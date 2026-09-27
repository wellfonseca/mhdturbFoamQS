# Como publicar este repositório (GitHub/GitLab + Zenodo)

Este repositório já está inicializado localmente com o primeiro commit. Falta
criá-lo num serviço público e arquivá-lo com DOI. Nada disto pode ser feito
automaticamente daqui — requer as suas contas.

## 1. Criar o repositório remoto

**GitHub** (via `gh`, se tiver o CLI autenticado):

```bash
cd mhdturbFoamQS
gh repo create mhdturbFoamQS --public --source=. --push \
   --description "Quasi-static MHD solver for OpenFOAM with non-uniform applied magnetic fields"
```

**GitLab** (ou manualmente pela interface web, criando um projecto vazio):

```bash
cd mhdturbFoamQS
git remote add origin git@gitlab.com:<utilizador>/mhdturbFoamQS.git
git branch -M main
git push -u origin main
```

Depois de saber o URL, actualize:

- `CITATION.cff` → campo `repository-code`
- `README.md` → secção 8 (Citation) e secção 7 do artigo

## 2. Arquivar no Zenodo (o DOI)

**Opção A — integração automática (recomendada).**

1. Entre em <https://zenodo.org> com a conta GitHub/GitLab.
2. *Settings → GitHub* (ou GitLab) → ligue o interruptor do repositório
   `mhdturbFoamQS`.
3. No GitHub, crie uma *release* com a etiqueta `v0.1.0`:
   ```bash
   git tag -a v0.1.0 -m "mhdturbFoamQS 0.1.0 — verification against Hartmann and Shercliff–Gold solutions"
   git push origin v0.1.0
   ```
4. O Zenodo arquiva automaticamente essa release e emite um **DOI**. Cole-o no
   `README.md`, no `CITATION.cff` (`doi:`), no artigo e no artigo companheiro.

**Opção B — envio manual.** Crie um tarball e carregue-o em
<https://zenodo.org/deposit/new>:

```bash
git archive --format=zip --prefix=mhdturbFoamQS-0.1.0/ -o /tmp/mhdturbFoamQS-0.1.0.zip HEAD
```

Metadados a preencher: título e resumo (copiar de `CITATION.cff`), tipo
*Software*, licença **GPL-3.0-or-later**, e os autores. Nos *related works*,
acrescente as duas obras de que deriva (FOSSEE/Radhakrishnan 2019 e
Tassone 2016) — está tudo no `NOTICE`.

## 3. Lista de verificação antes de submeter o artigo

- [ ] **A validação correu e passa.** `cd cases/channelHartmann && bash Allrun`
      tem de reproduzir a Tabela 2 do artigo com desvio < 1 % em todos os `Ha`.
      Enquanto isto não acontecer, o artigo não é submetível.
- [ ] Tabela 2 (canal plano) preenchida.
- [ ] Tabela 3 (tubo periódico vs Gold/Shercliff) preenchida.
- [ ] Tabelas 4 e 5 (demonstração com ímanes alternados) preenchidas.
- [ ] Secção 5.5 (a armadilha da parede condutora) demonstrada com números.
- [ ] URL do repositório e DOI do Zenodo inseridos no artigo.
- [ ] **Atribuições** confirmadas nos dois artigos: OpenFOAM, FOSSEE /
      R. Radhakrishnan (2019) e Tassone (2016).
- [ ] Texto obrigatório da CAPES com o ROR `00x0ma614` na secção *Funding*
      (a APC da Meccanica é custeada pela CAPES — acordo Springer Nature,
      revista híbrida confirmada na lista oficial).
- [ ] ORCID registado em <https://meusdados.capes.gov.br/>.
- [ ] Licença CC BY escolhida na submissão (exigência da CAPES).
- [ ] Carta de apresentação declara a existência do artigo companheiro
      (o das simulações) e a separação de conteúdos.

## 4. Se quiser publicar sem custo de APC

A *OpenFOAM Journal* publica exactamente este tipo de desenvolvimento **sem
taxa**, mas não tem factor de impacto. Mantenha-a como plano B caso a Meccanica
recuse.
