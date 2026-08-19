"""
Call Qwen VL model to extract structured information from paper figures.

Usage:
  python3 figure_extract.py <image_path> --mode genes --context ctx.txt
  python3 figure_extract.py <image_path> --prompt "what to extract"
  python3 figure_extract.py <image_path> --mode cell_counts
  python3 figure_extract.py <image_path> --mode json --prompt "..."
  echo "prompt" | python3 figure_extract.py <image_path>

Defaults: model=qwen3.6-plus (82% GT accuracy), thinking=on, thinking_budget=4000.
Use --no-thinking to disable thinking (qwen-vl-max mode).
The script reads Qwen_API_KEY from the project .env file by default.
Set DASHSCOPE_API_KEY env var to override.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import dashscope


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent


def load_env():
    env = {}
    env_path = PROJECT_ROOT / ".env"
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            env[key.strip()] = value.strip()
    return env


def get_api_key():
    key = os.getenv("DASHSCOPE_API_KEY")
    if key:
        return key
    env = load_env()
    return env.get("Qwen_API_KEY", "")


# Pre-built prompts for common extraction tasks
# Key design principles (from Qwen-VL experiments):
#   - OCR-only prompts suppress hallucination; structured categorization triggers "completion" behavior
#   - Paper context + figure-specific UP/DOWN rules dramatically improve accuracy
#   - qwen3.6-plus + thinking (82% GT) >> qwen-vl-max (43% GT)
#   - Flat gene→direction output is universal across papers; only Paper Context + How To Determine
#     sections are paper-specific (filled via {context} placeholder)
PROMPTS = {
    "genes": (
        "{context}\n\n"
        "TASK: Read ALL visible text from this dotplot image. For each gene, also record its "
        "expression trend (Up or Down) after using drug. Do NOT interpret or organize — just "
        "report what you literally see.\n\n"
        "Output this exact JSON:\n"
        "{\n"
        '  {"GENE1": "UP" or "DOWN"},\n'
        '  {"GENE2": "UP" or "DOWN"},\n'
        "  ...\n"
        "}\n\n"
        "Rules:\n"
        "- Report gene symbols left to right, exactly as printed. If a symbol has ambiguous "
        "characters, mark it like \"PAX[3/8?]\" but still include it.\n"
        "- Do NOT skip partially readable text — include it with uncertainty notes.\n"
        "- Do NOT add genes you cannot see. Better to miss 1 gene than hallucinate 10."
    ),
    "cell_counts": (
        "Extract all per-condition cell counts (n=...) visible in this figure. "
        "Look for counts in axis labels, figure captions embedded in the image, or panel titles. "
        "Output format: JSON object mapping condition name to cell count, "
        "e.g. {\"C\": 1234, \"T1\": 567, \"T2\": 432}. "
        "Also report the total if shown. Only report numbers you can clearly read."
    ),
    "heatmap_genes": (
        "This is a heatmap or dotplot from a scRNA-seq paper. "
        "For each row (gene) or column cluster that is labeled, list the gene names you can read. "
        "Pay special attention to: row labels on the y-axis, gene names next to dots/points, "
        "and any cluster/module labels (e.g. 'Module A', 'Cluster 1'). "
        "Output as JSON: {\"rows\": [\"gene1\",\"gene2\"], \"clusters\": {\"Module A\": [\"gene3\",\"gene4\"]}}. "
        "Only report text you can clearly read — do not guess partially visible labels."
    ),
    "violin_cell_counts": (
        "This is a violin plot or similar distribution plot. "
        "Extract the sample/condition names from the x-axis and any n=... cell counts "
        "visible above each violin or in the axis labels. "
        "Output as JSON: {\"condition\": n_cells, ...}. "
        "Only report numbers you can clearly read."
    ),
    "figure_overview": (
        "Describe this figure in detail. What type of plot is it? "
        "What are the x and y axes? What conditions or samples are shown? "
        "What genes, scores, or metrics are displayed? "
        "What is the main biological message? "
        "Be specific — quote labels and values you can read."
    ),
}


def call_qwen_vl(image_path: str, prompt: str, model: str = "qwen3.6-plus", api_key: str = None, enable_thinking: bool = True, thinking_budget: int = 4000):
    """Call Qwen VL model to analyze an image. Returns the model's text response."""
    if api_key is None:
        api_key = get_api_key()

    if not api_key:
        raise RuntimeError("No API key found. Set Qwen_API_KEY in .env or DASHSCOPE_API_KEY env var.")

    abs_path = Path(image_path).resolve()
    if not abs_path.exists():
        raise FileNotFoundError(f"Image not found: {image_path}")

    image_uri = f"file://{abs_path}"

    messages = [
        {
            "role": "user",
            "content": [
                {"image": image_uri},
                {"text": prompt},
            ],
        }
    ]

    extra_kwargs = {}
    if enable_thinking:
        extra_kwargs["enable_thinking"] = True
        if thinking_budget > 0:
            extra_kwargs["thinking_budget"] = thinking_budget

    resp = dashscope.MultiModalConversation.call(
        api_key=api_key,
        model=model,
        messages=messages,
        **extra_kwargs,
    )

    if resp.status_code != 200:
        raise RuntimeError(f"API error: code={resp.code}, message={resp.message}")

    return resp.output.choices[0].message.content[0]["text"]


def main():
    parser = argparse.ArgumentParser(
        description="Extract information from paper figures using Qwen VL model"
    )
    parser.add_argument(
        "image", nargs="?", default=None,
        help="Path to the figure image file"
    )
    parser.add_argument(
        "-p", "--prompt", default=None,
        help="Custom prompt describing what to extract from the image"
    )
    parser.add_argument(
        "-m", "--model", default="qwen3.6-plus",
        help="Qwen VL model name (default: qwen3.6-plus; qwen-vl-max for no-thinking mode)"
    )
    parser.add_argument(
        "--mode", default=None,
        choices=list(PROMPTS.keys()),
        help="Pre-built extraction mode (overrides --prompt)"
    )
    parser.add_argument(
        "--api-key", default=None,
        help="DashScope API key (default: read from .env)"
    )
    parser.add_argument(
        "--json", action="store_true",
        help="Output raw JSON with full response metadata"
    )
    parser.add_argument(
        "--list-modes", action="store_true",
        help="List available pre-built extraction modes and exit"
    )
    parser.add_argument(
        "--thinking", action="store_true", default=True,
        help="Enable deep thinking mode (default: on; use --no-thinking to disable)"
    )
    parser.add_argument(
        "--no-thinking", action="store_false", dest="thinking",
        help="Disable deep thinking mode"
    )
    parser.add_argument(
        "--thinking-budget", type=int, default=4000,
        help="Max thinking tokens when --thinking is enabled (default: 4000)"
    )
    parser.add_argument(
        "--context", default=None,
        help="Paper-specific context for the prompt (file path or inline string). "
             "Supports {context} placeholder in pre-built prompts."
    )

    args = parser.parse_args()

    if args.list_modes:
        print("Available modes:")
        for name, desc in PROMPTS.items():
            print(f"  {name}: {desc[:100]}...")
        return

    if not args.image:
        parser.print_help()
        sys.exit(1)

    # Determine prompt
    if args.mode:
        prompt = PROMPTS[args.mode]
    elif args.prompt:
        prompt = args.prompt
    elif not sys.stdin.isatty():
        prompt = sys.stdin.read().strip()
    else:
        prompt = PROMPTS["figure_overview"]

    # Handle --context: fill {context} placeholder or prepend
    context_text = ""
    if args.context:
        ctx_path = Path(args.context)
        if ctx_path.exists():
            context_text = ctx_path.read_text().strip()
        else:
            context_text = args.context.strip()

    if context_text:
        if "{context}" in prompt:
            prompt = prompt.replace("{context}", context_text)
        else:
            prompt = context_text + "\n\n" + prompt

    try:
        result = call_qwen_vl(
            image_path=args.image,
            prompt=prompt,
            model=args.model,
            api_key=args.api_key,
            enable_thinking=args.thinking,
            thinking_budget=args.thinking_budget,
        )

        if args.json:
            print(json.dumps({"status": "ok", "model": args.model, "result": result},
                             ensure_ascii=False, indent=2))
        else:
            print(result)

    except Exception as e:
        if args.json:
            print(json.dumps({"status": "error", "error": str(e)}, ensure_ascii=False))
        else:
            print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
