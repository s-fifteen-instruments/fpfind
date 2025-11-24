| Dataset | Description |
|---------|-------------|
| timestamp7_c14_rollover.Aa1fX.ts | 1000 events, channel 4 delayed copy of channel 1, 0b1100 filler, 862 actual events (138 dummy) |

Internal test data is referenced as a git submodule. To retrieve the submodule:

```bash
git submodule update --recursive --init    # first retrieval
git submodule update --recursive --remote  # subsequent retrievals
```
