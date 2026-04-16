from pathlib import Path
import csv
import re
from collections import defaultdict

MODULES_DIR = Path("modules")
TEMPLATE_PATH = Path("docs/README_template.md")
METADATA_PATH = Path("docs/module_seq_types.tsv")


def prompt(field_name: str) -> str:
    return input(f"{field_name}: ").strip()


def parse_version(version_str: str):
    return tuple(int(x) for x in re.findall(r"\d+", version_str))


def get_latest_version(module_path: Path) -> str:
    versions = [
        d.name for d in module_path.iterdir()
        if d.is_dir() and re.match(r"^\d+(\.\d+)*$", d.name)
    ]
    if not versions:
        return "UNKNOWN"
    return sorted(versions, key=parse_version)[-1]


def load_metadata():
    data = defaultdict(lambda: {
        "seq_types": set(),
        "purpose": None,
        "level": None,
        "input_format": None,
        "output_format": None,
    })

    with open(METADATA_PATH, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            module = row["module"]

            data[module]["seq_types"].add(row["seq_type"])
            data[module]["purpose"] = row["purpose"]
            data[module]["level"] = row["level"]
            data[module]["input_format"] = row["input_type"]
            data[module]["output_format"] = row["output_type"]

    for module in data:
        data[module]["seq_types"] = ", ".join(sorted(data[module]["seq_types"]))

    return data


def create_readme(module_path: Path, template: str, metadata: dict):
    readme_path = module_path / "README.md"

    if readme_path.exists():
        print(f"Skipping {module_path.name} (README exists)")
        return

    module_name = module_path.name
    version = get_latest_version(module_path)

    print(f"\n--- {module_name} (v{version}) ---")

    meta = metadata.get(module_name)

    if not meta:
        print("No metadata found — falling back to manual input")
        level = prompt("Level")
        module_purpose = prompt("Module purpose")
        input_format = prompt("Input format")
        output_format = prompt("Output format")
        seq_types = prompt("Seq types (comma separated)")
    else:
        level = meta["level"]
        module_purpose = meta["purpose"]
        input_format = meta["input_format"]
        output_format = meta["output_format"]
        seq_types = meta["seq_types"]

        print(f"(auto) level={level}, input={input_format}, output={output_format}, seq_types={seq_types}")

    # Only thing left to ask
    tool_name = prompt("Tool name")

    content = template.format(
        module_name=module_name,
        version=version,
        level=level,
        input_format=input_format,
        tool_name=tool_name,
        module_purpose=module_purpose,
        output_format=output_format,
        seq_types=seq_types,
    )

    readme_path.write_text(content)
    print(f"Created {readme_path}")

def main():
    if not MODULES_DIR.exists():
        print("modules/ directory not found")
        return

    template = TEMPLATE_PATH.read_text()
    metadata = load_metadata()

    for item in sorted(MODULES_DIR.iterdir()):
        if item.is_dir():
            create_readme(item, template, metadata)


if __name__ == "__main__":
    main()
