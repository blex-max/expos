# Extras

## **estimate-entropy**

This repo also contains a daughter tool for estimating the entropy rate of strings
in the same way as is done for assessing reference complexity. Useful if you're interested
in assessing the complexity of sequence data from any source (or any string at all!).

Enable compilation via:

```bash
cmake .. -DMAKE_DAUGHTER=ON
```

Then proceed with compilation as normal. See `./estimate-entropy --help` for usage.
