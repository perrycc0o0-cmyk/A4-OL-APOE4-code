# GitHub release procedure

Create the remote repository without adding a README, `.gitignore`, or license,
because the local repository already contains the curated files.

From the repository root:

```bash
git init -b main
git add .
git commit -m "Initial public analysis-code release"
git remote add origin https://github.com/OWNER/A4-OL-APOE4-code.git
git push -u origin main
```

Before making the repository public:

1. complete `docs/release_checklist.md`;
2. confirm the repository contains no private human metadata or reviewer tokens;
3. retain the approved MIT `LICENSE` file;
4. add the final `CITATION.cff`;
5. tag the exact manuscript version, for example `v1.0.0`;
6. create a GitHub release from the tag; and
7. archive that release in Zenodo or an equivalent repository and record the DOI.

For double-anonymized review, follow the journal's current anonymous-code review
instructions rather than exposing author identities through the public GitHub
profile during review.
