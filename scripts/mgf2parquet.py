import logging
import os
import pandas as pd

# Configure logger (accessible by caller)
logger = logging.getLogger(__name__)

def mgf_to_parquet(param):
    if os.path.exists(param['output_path']):
        logger.info(f"{param['output_path']} already exists")
        return 1
    if not os.path.exists(param['mgf_path']):
        logger.error(f"Missing {param['mgf_path']}")
        return 2
    try:
        records = []
        current = None
        field_config = param['field_type_dict']
        with open(param['mgf_path'], "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("COM="):
                    # Global comment, ignore
                    continue

                if line == "BEGIN IONS":
                    current = {"mz": [], "intensity": []}
                    continue

                if line == "END IONS":
                    if current:
                        # Ensure mz/intensity are at least empty lists
                        for k in ["mz", "intensity"]:
                            current[k] = current[k] if current[k] else []
                        records.append(current)
                        current = None
                    continue

                if current is None:
                    continue

                if "=" in line:  # Metadata line
                    key, value = line.split("=", 1)
                    key = key.strip().upper()
                    if key in field_config:
                        if key == "CHARGE" and field_config[key] == 'int':
                            val = value.strip()
                            try:
                                if val.endswith("+") or val.endswith("-"):
                                    current[key] = int(val[:-1]) * (1 if val.endswith("+") else -1)
                                else:
                                    current[key] = int(val)  # Handle unsigned case as well
                            except ValueError:
                                current[key] = val
                        # Convert to the specified data type
                        else:
                            dtype = field_config[key]
                            if dtype == "int":
                                try:
                                    current[key] = int(value)
                                except ValueError:
                                    pass
                            elif dtype == "float":
                                try:
                                    current[key] = float(value.split()[0])  # e.g., "380.18 12345"
                                except ValueError:
                                    pass
                            else:  # Default to string
                                current[key] = value.strip()
                else:  # mz-intensity line
                    try:
                        mz, intensity = map(float, line.split())
                        current["mz"].append(mz)
                        current["intensity"].append(intensity)
                    except ValueError:
                        pass  # Skip malformed lines

        # Build DataFrame
        df = pd.DataFrame(records)
        # Convert to Parquet
        df.to_parquet(param['output_path'])
        logger.info(f"{param['output_path']} has been successfully generated")
        return 0
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return -1
