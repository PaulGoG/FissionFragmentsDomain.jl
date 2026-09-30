# Dormant automation

The CI workflow and the Dependabot configuration are held here, outside the paths GitHub
reads, until the repository is public. Restoring them is two moves:

```
git mv .github/dormant/CI.yml .github/workflows/CI.yml
git mv .github/dormant/dependabot.yml .github/dependabot.yml
```

The gate the workflow applies is `julia check.jl --check`, which runs locally at any time.
