"""Guard the Claude Code integration: frontmatter must parse and reference real files.

Broken frontmatter fails silently at session start, so it is worth a test.
"""
import json
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")
ROOT = Path(__file__).resolve().parents[1]
AGENTS = sorted((ROOT / ".claude/agents").glob("*.md"))
COMMANDS = sorted((ROOT / ".claude/commands").glob("*.md"))


def _frontmatter(p: Path) -> dict:
    text = p.read_text()
    assert text.startswith("---"), f"{p.name} has no frontmatter"
    return yaml.safe_load(text.split("---")[1])


@pytest.mark.parametrize("p", AGENTS, ids=lambda p: p.name)
def test_agent_frontmatter(p):
    fm = _frontmatter(p)
    assert fm.get("name"), "name is required"
    assert fm.get("description"), "description drives delegation and is required"
    assert fm["name"] == p.stem, "name must match filename for @-mention to work"


@pytest.mark.parametrize("p", COMMANDS, ids=lambda p: p.name)
def test_command_frontmatter(p):
    fm = _frontmatter(p)
    assert fm.get("description"), "without description the picker shows the first body line"


def test_agents_exist():
    names = {p.stem for p in AGENTS}
    assert {"function-seed-researcher", "microbiologist-critic", "pipeline-runner"} <= names


def test_critic_is_read_only():
    """The critic must not be able to edit the seed it reviews."""
    fm = _frontmatter(ROOT / ".claude/agents/microbiologist-critic.md")
    tools = {t.strip() for t in str(fm.get("tools", "")).split(",")}
    assert not ({"Write", "Edit"} & tools), "critic must stay read-only"


def test_referenced_skill_exists():
    fm = _frontmatter(ROOT / ".claude/agents/function-seed-researcher.md")
    skill = fm.get("skills")
    assert skill, "seed researcher should declare its skill"
    name = skill if isinstance(skill, str) else skill[0]
    assert (ROOT / ".claude/skills" / name / "SKILL.md").exists()


def test_skill_scripts_resolve_repo_root():
    """Scripts moved under .claude/ must still find src/ regardless of depth."""
    for s in (ROOT / ".claude/skills/function-gene-seed-agent/scripts").glob("*.py"):
        assert "_repo_root" in s.read_text(), f"{s.name} uses a brittle relative path"


def test_settings_json_valid():
    s = json.loads((ROOT / ".claude/settings.json").read_text())
    assert "permissions" in s
