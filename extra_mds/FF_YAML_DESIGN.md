# Force Field YAML Format Design Document

## Problem Analysis

### Current Issues:
1. **YAML Builder** generates extended SMARTS: `[C;D4;H3][C;D4;H2]` (two atoms)
2. **Calculator** only types the FIRST atom in multi-atom patterns
3. **Result**: Atoms don't get typed → fallback to `generic_C`, `generic_O`

### Root Cause:
Inconsistency between pattern generation and pattern matching strategy.

## Ideal Solution: BASE SMARTS (Single-Atom Patterns)

### Rationale:
1. **Simplicity**: Each pattern describes ONE atom's local environment
2. **Completeness**: Every atom can be matched independently
3. **Efficiency**: No need to check connectivity during typing
4. **Standard**: Matches how OPLS-AA and other force fields work

### Format Specification:

```yaml
atom_types:
  - smarts: '[C;D4;H3]'      # Carbon with 4 neighbors, 3 of them H
    type_name: C_CH3         # Methyl carbon
    charge: 0.0
    sigma: 0.35
    epsilon: 0.276
    
  - smarts: '[C;D4;H2]'      # Carbon with 4 neighbors, 2 of them H  
    type_name: C_CH2         # Methylene carbon
    charge: 0.0
    sigma: 0.35
    epsilon: 0.276
    
  - smarts: '[O;D2;H1]'      # Oxygen with 2 neighbors, 1 of them H
    type_name: O_OH          # Hydroxyl oxygen
    charge: 0.0
    sigma: 0.296
    epsilon: 0.879
    
  - smarts: '[#1;D1]'        # Hydrogen with 1 neighbor (any H)
    type_name: H_all         # Generic hydrogen
    charge: 0.0
    sigma: 0.25
    epsilon: 0.126
```

### SMARTS Pattern Components:
- `[element]` - Element symbol or `#N` for atomic number
- `;DN` - Degree (number of bonded neighbors)
- `;HN` - Number of hydrogen atoms
- `;X1` - Connectivity (total bonds including H)
- `;+1` / `;-1` - Formal charge

### Examples:
```
[C;D4;H3]  → CH3 methyl carbon (sp3 with 3 H)
[C;D4;H2]  → CH2 methylene carbon (sp3 with 2 H)
[C;D4;H1]  → CH  methine carbon (sp3 with 1 H)
[C;D4;H0]  → C   quaternary carbon (sp3 with 0 H)
[C;D3;H1]  → CH  aromatic or sp2 carbon with 1 H
[C;D2;H0]  → C   aromatic carbon with no H
[#1;D1]    → H   any hydrogen
[O;D2;H1]  → OH  hydroxyl oxygen
[O;D2;H0]  → O   ether oxygen
[N;D3;H2]  → NH2 primary amine nitrogen
```

## Implementation Plan:

### Phase 1: Update YAML Builder (app.py)
1. ✅ Already groups by `base_smarts`
2. ✅ Already stores `base_smarts` in output
3. ❌ **FIX NEEDED**: Ensure 'smarts' field contains BASE pattern, not extended

### Phase 2: Update Calculator (calculator.py)  
1. ✅ Already reads 'smarts' field from YAML
2. ✅ Already matches with `GetSubstructMatches()`
3. ❌ **FIX NEEDED**: Handle single-atom patterns correctly (match[0] should work)
4. ❌ **VERIFY**: Ensure all atoms get typed before fallback

### Phase 3: Testing
1. Re-generate ethanol.yaml with base SMARTS
2. Verify all 9 atoms get typed correctly
3. No generic fallbacks should appear
4. Energy calculation should work

## Expected Output for Ethanol:

```yaml
atom_types:
  - smarts: '[C;D4;H3]'
    type_name: C_CH3
    # ... parameters
    
  - smarts: '[C;D4;H2]'
    type_name: C_CH2
    # ... parameters
    
  - smarts: '[O;D2;H1]'
    type_name: O_OH
    # ... parameters
    
  - smarts: '[#1;D1]'
    type_name: H_all
    # ... parameters
```

Result: 4 unique atom types covering all 9 atoms
- 1 methyl C (CH3)
- 1 methylene C (CH2)
- 1 hydroxyl O (OH)
- 6 hydrogens (all identical)
