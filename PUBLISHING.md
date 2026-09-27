# Publishing this repository (GitHub/GitLab + Zenodo)

The repository is already initialised locally and pushed. What remains is the
permanent archive with a DOI, which the target journals require.

## 1. Repository

The project lives at <https://github.com/wellfonseca/mhdturbFoamQS>.

To create the remote from scratch elsewhere (GitLab, or a fork):

```bash
gh repo create mhdturbFoamQS --public --source=. --push \
   --description "Quasi-static MHD solver for OpenFOAM with non-uniform applied magnetic fields"
```

or, for GitLab:

```bash
git remote add origin git@gitlab.com:<user>/mhdturbFoamQS.git
git branch -M main
git push -u origin main
```

Remember to update `repository-code` in `CITATION.cff` afterwards.

## 2. Archiving in Zenodo (the DOI)

**Option A — automatic integration (recommended).**

1. Sign in to <https://zenodo.org> **with GitHub**.
2. Go to <https://zenodo.org/account/settings/github/> and flip the switch for
   `wellfonseca/mhdturbFoamQS`.
3. On GitHub, publish a release for the existing tag `v0.1.0`:
   <https://github.com/wellfonseca/mhdturbFoamQS/releases/new?tag=v0.1.0>
4. Zenodo archives the release automatically and mints a DOI. The deposit
   metadata is pre-filled from `.zenodo.json`; paste the resulting DOI into
   `CITATION.cff`, `README.md` and the accompanying paper.

**Option B — manual upload.**

Create an archive from a tag and upload it at <https://zenodo.org/deposit/new>:

```bash
git archive --format=zip --prefix=mhdturbFoamQS-0.1.0/ \
    -o mhdturbFoamQS-0.1.0.zip v0.1.0
```

Metadata: copy from `.zenodo.json` (title, description, creators, keywords),
resource type *Software*, licence **GPL-3.0-or-later**, and add the two upstream
works as related identifiers — see `NOTICE`.

## 3. Checklist before submitting the paper

- [ ] **Verification passes.** `cd cases/channelHartmann && bash Allrun` must
      reproduce Table 2 of the paper with a deviation below 1 % at every `Ha`.
      Until this holds, the paper cannot be submitted.
- [ ] Tables 2 and 3 (plane channel and circular pipe) filled in.
- [ ] Tables 4 and 5 (alternating-magnet demonstration) filled in.
- [ ] Section 5.4 (the conducting-wall trap) demonstrated numerically.
- [ ] Repository URL and Zenodo DOI inserted in the paper.
- [ ] **Attribution** confirmed in both papers: OpenFOAM, FOSSEE /
      R. Radhakrishnan (2019) and Tassone (2016).
- [ ] Mandatory CAPES funding statement, with ROR `00x0ma614` (the Meccanica
      article processing charge is covered by CAPES under the Springer Nature
      agreement — a hybrid journal, confirmed in the official list).
- [ ] ORCID registered at <https://meusdados.capes.gov.br/>.
- [ ] CC BY licence selected at submission (CAPES requirement).
- [ ] Cover letter states the existence of the companion paper (the simulation
      study) and the separation of contents.

## 4. If no article processing charge can be paid

The *OpenFOAM Journal* publishes exactly this kind of development **free of
charge**, but has no impact factor. Keep it as a fallback should Meccanica
decline.
